set no fetch;

store(
  redimension(
    substitute(
      equi_join(
        cast(
          $NAMESPACE.ASSOC_${ASSOCIATION}_FIELD,
          (description, description_old),
          (notes,       notes_old)),
        project(
          apply(
            aio_input(
              '/opt/biobankengine/GlobalBioBankEngineRepo/gbe_data/$FILE',
              num_attributes:6,
              header:1),
            title_match,            a1,
            description,            a3,
            notes,                  'number_of_participants=' + a2 +
                                    ';short_name=' + a4 +
                                    ';category=' + a0),
          title_match,
          description,
          notes),
        left_names:title,
        right_names:title_match,
        keep_dimensions:true,
        left_outer:true),
      build(<subval:string not null>[i=0:0], 'Uploaded on $DATE'),
      notes),
    <title                  : string,
     value_type             : string,
     description            : string,
     notes                  : string NOT NULL,
     encoding_id            : int64,
     pvalue_threshold       : double>
    [field_id=0:*:0:100000]),
  $NAMESPACE.ASSOC_${ASSOCIATION}_FIELD_UPDATED);
