set no fetch;

insert(
  redimension(
    project(
      apply(
        aio_input(
          '/opt/biobankengine/GlobalBioBankEngineRepo/gbe_data/bbj_variant_annots_gbe.scidb.$I',
          num_attributes:31),
        variant_id,    int64(a0),
        chrom,         int64(a1),
        pos,           int64(a2),
        ref,           a3,
        alt,           a4,
        rsid,          iif(substr(a5, 0, 2) = 'rs', a5, null),
        id,            iif(substr(a5, 0, 2) = 'rs',
                           dcast(substr(a5, 2, strlen(a5)), int64(null)),
                       null),
        gene,          a6,
        consequence,   a7,
        impact,        iif(a13 = '""', null, rsub(a13, 's/([^|]*\|){2}([^|]*)\|.*/$2/')),
        hgvsp,         iif(a8 = '""', null, a8),
        lof,           iif(a9 = '""', null, a9),
        annotations,   iif(a13 = '""', null, a13),
        miss,          dcast(a14, double(null)),
        miss_bileve,   dcast(a15, double(null)),
        miss_wcsg,     dcast(a16, double(null)),
        hwep,          dcast(a17, double(null)),
        maf,           dcast(a18, double(null)),
        ld,            iif(a19 = 'False', false,
                           iif(a19 = 'True', true, null)),
        wcsg_only,     iif(a20 = 'False', false,
                           iif(a20 = 'True', true, null)),
        bileve_only,   iif(a21 = 'False', false,
                           iif(a21 = 'True', true, null)),
        missingness,   iif(a23 = '0', false,
                           iif(a23 = '1', true, null)),
        hwe,           iif(a24 = '0', false,
                           iif(a24 = '1', true, null)),
        mcpi,          iif(a25 = '0', false,
                           iif(a25 = '1', true, null)),
        gnomad_filter, iif(a26 = '""', null, a26),
        mgi,           iif(a27 = '""', null, a27),
        mgi_notes,     iif(a28 = '""', null, a28),
        all_filters,   dcast(a29, int8(null))
        -- gene_symbol,   iif(a30 = '' or a30 = '""', null,
        --                    iif(substr(a30, 0, 2) = '""',
        --                        substr(a30, 2, strlen(a30)), a30))
      ),
      variant_id,
      chrom,
      pos,
      ref,
      alt,
      rsid,
      id,
      gene,
      consequence,
      impact,
      hgvsp,
      lof,
      annotations,
      miss,
      miss_bileve,
      miss_wcsg,
      hwep,
      maf,
      ld,
      wcsg_only,
      bileve_only,
      missingness,
      hwe,
      mcpi,
      gnomad_filter,
      mgi,
      mgi_notes,
      all_filters
      -- gene_symbol
    ),
    LOAD.VARIANT),
  LOAD.VARIANT);
