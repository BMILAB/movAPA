library���ƺ�����genome(BSgenome��ʽ)�Ķ�Ӧ��ϵ
library	BSgenomeObjname
BSgenome.Athaliana.ENSEMBL.TAIR10	Arabidopsis
BSgenome.Oryza.ENSEMBL.IRGSP1	Oryza
BSgenome.Indica.ENSEMBL.ASM465V1	Indica
BSgenome.Medicago.ENSEMBL.MEDTRA17	Medicago
BSgenome.Trifolium.ENSEMBL.TRPR	Trifolium
BSgenome.Chlamydomonas.ENSEMBL.V5	Chlamydomonas
BSgenome.Bamboo.NCGI.V1	Bamboo
##�����
##7�����ִ����BSgeome��chr�ĳ��Ⱥ�ԭʼfasta����һ��


���ؽ���Щ��install���������library: �½�һ��new project > existing directory ѡ��һ������ > �Ҳ�build...build and ... �����ɰ�
��ATH���ذ�װ���ϱ�ġ�BSgenomeObjname������ͬ����genome����
library(BSgenome.Athaliana.ENSEMBL.TAIR10)
seqlengths(Arabidopsis)
getSeq(Arabidopsis,"1",start=1,end=100,strand="+")
