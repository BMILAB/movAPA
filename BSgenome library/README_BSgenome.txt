library名称和物种genome(BSgenome格式)的对应关系
library	BSgenomeObjname
BSgenome.Athaliana.ENSEMBL.TAIR10	Arabidopsis
BSgenome.Oryza.ENSEMBL.IRGSP1	Oryza
BSgenome.Indica.ENSEMBL.ASM465V1	Indica
BSgenome.Medicago.ENSEMBL.MEDTRA17	Medicago
BSgenome.Trifolium.ENSEMBL.TRPR	Trifolium
BSgenome.Chlamydomonas.ENSEMBL.V5	Chlamydomonas
BSgenome.Bamboo.NCGI.V1	Bamboo
##检查结果
##7个物种打包成BSgeome的chr的长度和原始fasta长度一致


本地将这些包install便可以生成library: 新建一个new project > existing directory 选择一个物种 > 右侧build...build and ... 制作成包
如ATH本地安装后，上表的“BSgenomeObjname”代表不同物种genome名称
library(BSgenome.Athaliana.ENSEMBL.TAIR10)
seqlengths(Arabidopsis)
getSeq(Arabidopsis,"1",start=1,end=100,strand="+")
