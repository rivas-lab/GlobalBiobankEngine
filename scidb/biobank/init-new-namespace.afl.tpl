create_namespace($NAMESPACE);

store(
  apply(
    build(<description:string>[i=0:0], '$NAMESPACE'),
    number_of_participants, 0,
    number_of_variants,     0,
    assembly,               'hg38'),
  $NAMESPACE.NAMESPACE_DESCRIPTION);

create array $NAMESPACE.ASSOCIATION_SET
  <name               : string,
   variant_namespace  : string NOT NULL,
   number_of_fields   : int64,
   number_of_variants : int64,
   date_modified      : datetime,
   note               : string>
   [id = 0:*:0:1];
