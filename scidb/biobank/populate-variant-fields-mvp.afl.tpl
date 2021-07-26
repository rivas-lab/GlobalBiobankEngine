
store(
  cast(
    apply(
      project(
        filter(
          attributes($NAMESPACE.VARIANT),
          name <> 'id' and name <> 'alt_id'),
        name, type_id),
      output_optional, iif(name='annotations', true, false)),
    <name:string, type:string, output_optional:bool>[i]),
  $NAMESPACE.VARIANT_FIELDS);
