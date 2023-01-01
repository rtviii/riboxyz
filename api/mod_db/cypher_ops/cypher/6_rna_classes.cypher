UNWIND [ 
  "5SrRNA"  ,
  "5.8SrRNA",
  "12SrRNA" ,
  "16SrRNA" ,
  "21SrRNA" ,
  "23SrRNA" ,
  "25SrRNA" ,
  "28SrRNA" ,
  "35SrRNA" ,
  "mRNA"    ,
  "tRNA"    ]  as rnaclass create (n:NomenclatureClass {class_id:rnaclass})