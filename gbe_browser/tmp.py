vep_annotation = str(variant['annotations'])
info_field = dict([(x.split('=', 1)) if '=' in x else (x, x) for x in re.split(';(?=\w)', vep_annotation)])
consequence_array = info_field['CSQ'].split(',') if 'CSQ' in info_field else []
annotations = [dict(zip(utils.vep_field_names, x.split('|'))) for x in consequence_array]
coding_annotations = [ann for ann in annotations if 'Feature' in ann and ann['Feature'].startswith('ENST')]
vep_annotations = [ann for ann in coding_annotations]
variant['vep_annotations'] = vep_annotations
variant['filter'] =  'PASS'
#if variant['all_filters'] != 1 and variant['all_filters'] != 2 and variant['all_filters'] != 3 else 'FAIL'
            # variant['vep_annotations'] = order_vep_by_csq(variant['vep_annotations'])  # Adds major_consequence
           # ordered_csqs = [x['major_consequence'] for x in variant['vep_annotations']]
           # ordered_csqs = reduce(lambda x, y: ','.join([x, y]) if y not in x else x, ordered_csqs, '').split(',') # Close but not quite there
           # consequences = defaultdict(lambda: defaultdict(list))
           # for annotation in variant['vep_annotations']:
           #     annotation['HGVS'] = get_proper_hgvs(annotation)
           #     consequences[annotation['major_consequence']][annotation['Gene']].append(annotation)
variant['vep_annotations'] = utils.order_vep_by_csq(variant['vep_annotations'])  # Adds major_consequence
ordered_csqs = [x['major_consequence'] for x in variant['vep_annotations']]
ordered_csqs = reduce(lambda x, y: ','.join([x, y]) if y not in x else x, ordered_csqs, '').split(',') # Close but not quite there
consequences = defaultdict(lambda: defaultdict(list))
for annotation in variant['vep_annotations']:
    annotation['HGVS'] = utils.get_proper_hgvs(annotation)
    consequences[annotation['major_consequence']][annotation['Gene']].append(annotation)
assocset = str(db.list_association_sets()['name'][0])
icdstats = db.get_association_data(association_set=assocset, chromosome=chrom, position=pos, association_fields = ('pvalue', 'title', 'beta', 'odds_ratio', 'beta', 'se'))
icdstats.query('ref == @refl and alt == @altl', inplace = True)
icdstats = icdstats.to_dict('records')
indexes = []
seend = {}
phef = db.get_phenotype_fields(association_set=str(db.list_association_sets()['name'][0]))
for idx in range(0,len(icdstats)):
    # ICD10=T81/Name=Complications of procedures, not elsewhere classified/Chapter=T/L95OR=0.97/U95OR=2.04/OR=1.40/pvalue=0.0756463/l10pval=1.12/Case=1229
    item = icdstats[idx]
    icd10 = str(item['title'])
    item['Code'] = icd10
    if 'title' not in item or (numpy.isnan(item['odds_ratio']) and numpy.isnan(item['beta'])) or icd10 == "HC90":
        item['Name'] = 'NA'
        item['Group'] = 'NA'
        item['OR'] = 1
        item['LOR'] = 0
        item['L95OR'] = 1
        item['U95OR'] = 1
        item['pvalue'] = 1
        item['l10pval'] = 0
        item['Case'] = 'NA'
        item['SE'] = 0
        indexes.append(idx)
    else:
        item['Name'] = str(phef[phef['title'] == icd10]['notes'].squeeze()).split(';')[1].split('=')[1]
        item['Case'] = str(phef[phef['title'] == icd10]['notes'].squeeze()).split(';')[0].split('=')[1]
        item['Group'] = icd10 # default value
        groups = ['RH', 'FH', 'HC', 'cancer', 'ADD', 'INI_FC', 'QT_FC', 'BIN_FC', 'INI', 'MED', 'BIN', 'BRMRI', 'BROADBIN', 'BROADQT']
        for groupidx in range(0, len(groups)):
            if icd10.startswith(groups[groupidx]):
                item['Group'] = groups[groupidx]
                break
            if numpy.isnan(item['odds_ratio']):
                item['OR'] = 1
            else:
                item['OR'] = float(item['odds_ratio'])
            if numpy.isnan(item['beta']):
                item['LOR'] = numpy.log(float(item['OR']))
            else:
                item['LOR'] = float(item['beta'])
            item['L95OR'] = float(item['LOR']) - 1.96*float(item['se'])
            item['U95OR'] = float(item['LOR']) + 1.96*float(item['se'])
            item['pvalue'] = float(item['pvalue'])
            item['l10pval'] = -numpy.log10(float(item['pvalue']))
            item['SE'] = float(item['se'])
            if float(item['pvalue']) == 0:
                item['pvalue'] = numpy.finfo(float).eps
                item['pvalue'] = float(item['pvalue'])
                item['l10pval'] = 250
                # item['Case'] = icd10info[0]['Case']
            item['log10pvalue'] = float(item['l10pval'])
            se =  float(item['se'])
            if item['l10pval'] < 1 or se >= .5 or (se >= .08 and item['OR'] == item['LOR']) or int(item['Case']) <= 100  or item['Code'] == "HC67" or icd10 in seend:
                indexes.append(idx)
            seend[icd10] = icd10
    for index in sorted(indexes, reverse=True):
        del icdstats[index]
        print('Rendering variant: %s' % variant_str)
        return render_template(
            'variant.html',
            variant=variant,
            icdstats=icdstats,
            consequences=consequences,
            ordered_csqs=ordered_csqs,
            debug_message='',
            namespace=namespace,
            #plot_pval_data = [{'x': [1, 2, 3], 'y': [1,2, 3]}],
            pval_slider_max = variant_page_pval_slider_max(icdstats),
            plot_pval_data = variant_page_plot_pval_data(icdstats),
            plot_lor_data  = variant_page_plot_lor_data(icdstats),
        )
    except Exception as e:
        print('Failed on variant:', variant_str, '; Error=', traceback.format_exc())
        abort(404)

def variant_page_data_prep_sub(icdstats, sort_key='log10pvalue'):
    plot_d_raw = collections.defaultdict(list)
    keys = icdstats[0].keys()
    for key in keys:
        plot_d_raw[key] = np.array([x[key] for x in icdstats])
    plot_df = pandas.DataFrame(plot_d_raw).sort_values(
        by=['Group', sort_key], ascending=[True, False]
    )
    plot_d_dict = collections.defaultdict(collections.defaultdict)
    groups = sorted(set(plot_df['Group']))
    for group in groups:
        for key in keys:
            plot_d_dict[group][key] = list(plot_df[plot_df['Group'] == group][key])
    for group in groups:
        for key in ['OR', 'LOR', 'L95OR', 'U95OR', 'pvalue', 'SE', 'log10pvalue']:
            plot_d_dict[group][key] = [float(x) for x in plot_d_dict[group][key]]
    for group in groups:
        #error_bar = {'L95OR': -1, 'U95OR': 1}
        #for key in error_bar.keys():
        #    diff = np.array(plot_d_dict[group][key]) - np.array(plot_d_dict[group]['LOR'])
        #    plot_d_dict[group]['d{}'.format(key)] = [0 if np.isnan(x) else np.abs(x) for x in diff]
        plot_d_dict[group]['196SE'] = list( 1.96 * np.array(plot_d_dict[group]['SE']) )
    for group in groups:
        if group in set(['INI', 'INI_FC', 'QT_FC', 'BROADQT']):
            beta_or_lor = 'BETA'
            beta_or_lor_val = plot_d_dict[group]['LOR']
            beta_or_lor_l95 = np.array(plot_d_dict[group]['LOR']) - 1.96 * np.array(plot_d_dict[group]['SE'])
            beta_or_lor_u95 = np.array(plot_d_dict[group]['LOR']) + 1.96 * np.array(plot_d_dict[group]['SE'])
        else:
            beta_or_lor = 'OR'
            beta_or_lor_val = np.exp(np.array(plot_d_dict[group]['LOR']))
            beta_or_lor_l95 = np.exp(np.array(plot_d_dict[group]['LOR']) - 1.96 * np.array(plot_d_dict[group]['SE']))
            beta_or_lor_u95 = np.exp(np.array(plot_d_dict[group]['LOR']) + 1.96 * np.array(plot_d_dict[group]['SE']))
        group_len = len(plot_d_dict[group]['Code'])
        plot_d_dict[group]['text'] = [
            '{}. Case: {}, P-value: {:.3e}, {} = {:.5f} (95% [{:.5f}, {:.5f}]), SE = {:.5f}'.format(
                ''.join([c if c != '_' else ' ' for c in x[0]]), x[1], x[2], x[3], x[4], x[5], x[6], x[7]
            ) for x in zip(
                plot_d_dict[group]['Name'],
                plot_d_dict[group]['Case'],
                plot_d_dict[group]['pvalue'],
                [beta_or_lor] * group_len,
                beta_or_lor_val,
                beta_or_lor_l95,
                beta_or_lor_u95,
                plot_d_dict[group]['SE'],
                #plot_d_dict[group]['L95OR'],
                #plot_d_dict[group]['U95OR'],
            )
        ]
    return plot_d_dict


def variant_page_pval_slider_max(icdstats):
#    return 100
    plot_d_dict = variant_page_data_prep_sub(icdstats, sort_key='log10pvalue')
    return np.max([np.max(v['log10pvalue']) for v in plot_d_dict.values()]) 


def variant_page_plot_pval_data(icdstats):
    plot_d_dict = variant_page_data_prep_sub(icdstats, sort_key='log10pvalue')
    groups = plot_d_dict.keys()
    plot_d = [{
        'x':    plot_d_dict[group]['Code'],
        'y':    plot_d_dict[group]['l10pval'],
        'text': plot_d_dict[group]['text'],
        'name': group,
        'type': 'scatter',
        'mode': 'markers',
        'marker': {'size': 16, },
        'hoverinfo':'x+text', 
    } for group in groups]
 
    return plot_d


def variant_page_plot_lor_data(icdstats):
    plot_d_dict = variant_page_data_prep_sub(icdstats, sort_key='log10pvalue')
    groups = plot_d_dict.keys()
    plot_d = [{
        'x':    plot_d_dict[group]['Code'],
        'y':    plot_d_dict[group]['LOR'],
        'error_y': {
            'type'       : 'data',
            'symmetric'  : 'false',
            'array'      : plot_d_dict[group]['196SE'],
            'arrayminus' : plot_d_dict[group]['196SE'],
        },
        'text': plot_d_dict[group]['text'],
        'name': group,
        'type': 'scatter',
        'mode': 'markers',
        'marker': {'size': 16, },
        'hoverinfo':'x+text', 
    } for group in groups]
 
    return plot_d


@app.route('/<namespace>/coding/<icd_str>')
def icd_page(namespace,icd_str):
#    if not check_credentials():
#        return redirect(url_for('login'))
    db = get_db(namespace)
        # icdlabel = str(icd_str)
        # .strip('ADD').strip('INI').strip('BRMRI')
        # Try several cutoffs to see if there are too many variants to
        # render.  Arbitrary 10k max.
        # passing = False
    cutoff = None
    icd = None
    for p in [.0001]:
        assocset = str(db.list_association_sets()['name'][0])
        df = db.get_phenotype_fields(association_set=assocset, include_pvalue_threshold=True)
        field_identifier = int(df[df['title'] == icd_str]['field_id'])
        pvalthr = float(df[df['title'] == icd_str]['pvalue_threshold'])
        if np.isnan(pvalthr):
            pvalthr = p
        pvalthr = max(5e-7, pvalthr)
