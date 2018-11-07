##
# File:    SequenceExamples.py
# Date:    2-Dec-2009
# Updates:
#
# 10-Jan-2010 jdw  Random and selected mutations in sequence data added for testing.
# 12-Feb-2010 jdw  Adopt SequenceFeature()
# 20-Apr-2010 jdw  Port to module seqmodule.
#
##
"""
Selected examples of author coordinate and reference sequences from PDB entry 3IJE.

Sequences are returned as a tuple containing (aa3, index, comment placeholder) -

"""
__docformat__ = "restructuredtext en"
__author__    = "John Westbrook"
__email__     = "jwest@rcsb.rutgers.edu"
__license__   = "Creative Commons Attribution 3.0 Unported"
__version__   = "V0.07"


import sys, string, random, copy
from wwpdb.apps.seqmodule.util.SequenceReferenceData import SequenceReferenceData
from wwpdb.apps.seqmodule.util.SequenceLabel         import SequenceFeature

class SequenceExamples(object):
    def __init__(self,idCode="3IJE",verbose=False,log=sys.stderr):
        self.__verbose=verbose
        self.__lfh=log
        self.__srd=SequenceReferenceData(verbose=self.__verbose,log=self.__lfh)
        
        self.idCode=idCode

        self.aaList=["ALA","ARG","ASN","ASP","ASX","CYS","GLN",
                     "GLU","GLX","GLY","HIS","ILE","LEU","LYS",
                     "MET","PHE","PRO","SER","THR","TRP","TYR","VAL","UNK"]

        self.seqIdList=['A','B']

        #
        # DBREF  3IJE A    1   967  UNP    P06756   ITAV_HUMAN      31    997             
        # DBREF  3IJE B    1   695  UNP    P05106   ITB3_HUMAN      27    721             
        # >sp|P06756|ITAV_HUMAN Integrin alpha-V OS=Homo sapiens GN=ITGAV PE=1 SV=2
        #
        self.refSequenceA='''MAFPPRRRLRLGPRGLPLLLSGLLLPLCRAFNLDVDSPAEYSGPEGSYFGFAVDFFVPSA
        SSRMFLLVGAPKANTTQPGIVEGGQVLKCDWSSTRRCQPIEFDATGNRDYAKDDPLEFKS
        HQWFGASVRSKQDKILACAPLYHWRTEMKQEREPVGTCFLQDGTKTVEYAPCRSQDIDAD
        GQGFCQGGFSIDFTKADRVLLGGPGSFYWQGQLISDQVAEIVSKYDPNVYSIKYNNQLAT
        RTAQAIFDDSYLGYSVAVGDFNGDGIDDFVSGVPRAARTLGMVYIYDGKNMSSLYNFTGE
        QMAAYFGFSVAATDINGDDYADVFIGAPLFMDRGSDGKLQEVGQVSVSLQRASGDFQTTK
        LNGFEVFARFGSAIAPLGDLDQDGFNDIAIAAPYGGEDKKGIVYIFNGRSTGLNAVPSQI
        LEGQWAARSMPPSFGYSMKGATDIDKNGYPDLIVGAFGVDRAILYRARPVITVNAGLEVY
        PSILNQDNKTCSLPGTALKVSCFNVRFCLKADGKGVLPRKLNFQVELLLDKLKQKGAIRR
        ALFLYSRSPSHSKNMTISRGGLMQCEELIAYLRDESEFRDKLTPITIFMEYRLDYRTAAD
        TTGLQPILNQFTPANISRQAHILLDCGEDNVCKPKLEVSVDSDQKKIYIGDDNPLTLIVK
        AQNQGEGAYEAELIVSIPLQADFIGVVRNNEALARLSCAFKTENQTRQVVCDLGNPMKAG
        TQLLAGLRFSVHQQSEMDTSVKFDLQIQSSNLFDKVSPVVSHKVDLAVLAAVEIRGVSSP
        DHVFLPIPNWEHKENPETEEDVGPVVQHIYELRNNGPSSFSKAMLHLQWPYKYNNNTLLY
        ILHYDIDGPMNCTSDMEINPLRIKISSLQTTEKNDTVAGQGERDHLITKRDLALSEGDIH
        TLGCGVAQCLKIVCQVGRLDRGKSAILYVKSLLWTETFMNKENQNHSYSLKSSASFNVIE
        FPYKNLPIEDITNSTLVTTNVTWGIQPAPMPVPVWVIILAVLAGLLLLAVLVFVMYRMGF
        FKRVRPPQEEQEREQLQPHENGEGNSET'''
        #
        # >sp|P05106|ITB3_HUMAN Integrin beta-3 OS=Homo sapiens GN=ITGB3 PE=1 SV=2
        self.refSequenceB='''
        MRARPRPRPLWATVLALGALAGVGVGGPNICTTRGVSSCQQCLAVSPMCAWCSDEALPLG
        SPRCDLKENLLKDNCAPESIEFPVSEARVLEDRPLSDKGSGDSSQVTQVSPQRIALRLRP
        DDSKNFSIQVRQVEDYPVDIYYLMDLSYSMKDDLWSIQNLGTKLATQMRKLTSNLRIGFG
        AFVDKPVSPYMYISPPEALENPCYDMKTTCLPMFGYKHVLTLTDQVTRFNEEVKKQSVSR
        NRDAPEGGFDAIMQATVCDEKIGWRNDASHLLVFTTDAKTHIALDGRLAGIVQPNDGQCH
        VGSDNHYSASTTMDYPSLGLMTEKLSQKNINLIFAVTENVVNLYQNYSELIPGTTVGVLS
        MDSSNVLQLIVDAYGKIRSKVELEVRDLPEELSLSFNATCLNNEVIPGLKSCMGLKIGDT
        VSFSIEAKVRGCPQEKEKSFTIKPVGFKDSLIVQVTFDCDCACQAQAEPNSHRCNNGNGT
        FECGVCRCGPGWLGSQCECSEEDYRPSQQDECSPREGQPVCSQRGECLCGQCVCHSSDFG
        KITGKYCECDDFSCVRYKGEMCSGHGQCSCGDCLCDSDWTGYYCNCTTRTDTCMSSNGLL
        CSGRGKCECGSCVCIQPGSYGDTCEKCPTCPDACTFKKECVECKKFDRGALHDENTCNRY
        CRDEIESVKELKDTGKDAVNCTYKNEDDCVVRFQYYEDSSGKSILYVVEEPECPKGPDIL
        VVLLSVMGAILLIGLAALLIWKLLITIHDRKEFAKFEEERARAKWDTANNPLYKEATSTF
        TNITYRGT'''

        self.authSequenceA='''PHE ASN LEU ASP VAL ASP SER PRO ALA GLU TYR SER GLY          
        PRO GLU GLY SER TYR PHE GLY PHE ALA VAL ASP PHE PHE          
        VAL PRO SER ALA SER SER ARG MET PHE LEU LEU VAL GLY          
        ALA PRO LYS ALA ASN THR THR GLN PRO GLY ILE VAL GLU          
        GLY GLY GLN VAL LEU LYS CYS ASP TRP SER SER THR ARG          
        ARG CYS GLN PRO ILE GLU PHE ASP ALA THR GLY ASN ARG          
        ASP TYR ALA LYS ASP ASP PRO LEU GLU PHE LYS SER HIS          
        GLN TRP PHE GLY ALA SER VAL ARG SER LYS GLN ASP LYS          
        ILE LEU ALA CYS ALA PRO LEU TYR HIS TRP ARG THR GLU          
        MET LYS GLN GLU ARG GLU PRO VAL GLY THR CYS PHE LEU          
        GLN ASP GLY THR LYS THR VAL GLU TYR ALA PRO CYS ARG          
        SER GLN ASP ILE ASP ALA ASP GLY GLN GLY PHE CYS GLN          
        GLY GLY PHE SER ILE ASP PHE THR LYS ALA ASP ARG VAL          
        LEU LEU GLY GLY PRO GLY SER PHE TYR TRP GLN GLY GLN          
        LEU ILE SER ASP GLN VAL ALA GLU ILE VAL SER LYS TYR          
        ASP PRO ASN VAL TYR SER ILE LYS TYR ASN ASN GLN LEU          
        ALA THR ARG THR ALA GLN ALA ILE PHE ASP ASP SER TYR          
        LEU GLY TYR SER VAL ALA VAL GLY ASP PHE ASN GLY ASP          
        GLY ILE ASP ASP PHE VAL SER GLY VAL PRO ARG ALA ALA          
        ARG THR LEU GLY MET VAL TYR ILE TYR ASP GLY LYS ASN          
        MET SER SER LEU TYR ASN PHE THR GLY GLU GLN MET ALA          
        ALA TYR PHE GLY PHE SER VAL ALA ALA THR ASP ILE ASN          
        GLY ASP ASP TYR ALA ASP VAL PHE ILE GLY ALA PRO LEU          
        PHE MET ASP ARG GLY SER ASP GLY LYS LEU GLN GLU VAL          
        GLY GLN VAL SER VAL SER LEU GLN ARG ALA SER GLY ASP          
        PHE GLN THR THR LYS LEU ASN GLY PHE GLU VAL PHE ALA          
        ARG PHE GLY SER ALA ILE ALA PRO LEU GLY ASP LEU ASP          
        GLN ASP GLY PHE ASN ASP ILE ALA ILE ALA ALA PRO TYR          
        GLY GLY GLU ASP LYS LYS GLY ILE VAL TYR ILE PHE ASN          
        GLY ARG SER THR GLY LEU ASN ALA VAL PRO SER GLN ILE          
        LEU GLU GLY GLN TRP ALA ALA ARG SER MET PRO PRO SER          
        PHE GLY TYR SER MET LYS GLY ALA THR ASP ILE ASP LYS          
        ASN GLY TYR PRO ASP LEU ILE VAL GLY ALA PHE GLY VAL          
        ASP ARG ALA ILE LEU TYR ARG ALA ARG PRO VAL ILE THR          
        VAL ASN ALA GLY LEU GLU VAL TYR PRO SER ILE LEU ASN          
        GLN ASP ASN LYS THR CYS SER LEU PRO GLY THR ALA LEU          
        LYS VAL SER CYS PHE ASN VAL ARG PHE CYS LEU LYS ALA          
        ASP GLY LYS GLY VAL LEU PRO ARG LYS LEU ASN PHE GLN          
        VAL GLU LEU LEU LEU ASP LYS LEU LYS GLN LYS GLY ALA          
        ILE ARG ARG ALA LEU PHE LEU TYR SER ARG SER PRO SER          
        HIS SER LYS ASN MET THR ILE SER ARG GLY GLY LEU MET          
        GLN CYS GLU GLU LEU ILE ALA TYR LEU ARG ASP GLU SER          
        GLU PHE ARG ASP LYS LEU THR PRO ILE THR ILE PHE MET          
        GLU TYR ARG LEU ASP TYR ARG THR ALA ALA ASP THR THR          
        GLY LEU GLN PRO ILE LEU ASN GLN PHE THR PRO ALA ASN          
        ILE SER ARG GLN ALA HIS ILE LEU LEU ASP CYS GLY GLU          
        ASP ASN VAL CYS LYS PRO LYS LEU GLU VAL SER VAL ASP          
        SER ASP GLN LYS LYS ILE TYR ILE GLY ASP ASP ASN PRO          
        LEU THR LEU ILE VAL LYS ALA GLN ASN GLN GLY GLU GLY          
        ALA TYR GLU ALA GLU LEU ILE VAL SER ILE PRO LEU GLN          
        ALA ASP PHE ILE GLY VAL VAL ARG ASN ASN GLU ALA LEU          
        ALA ARG LEU SER CYS ALA PHE LYS THR GLU ASN GLN THR          
        ARG GLN VAL VAL CYS ASP LEU GLY ASN PRO MET LYS ALA          
        GLY THR GLN LEU LEU ALA GLY LEU ARG PHE SER VAL HIS          
        GLN GLN SER GLU MET ASP THR SER VAL LYS PHE ASP LEU          
        GLN ILE GLN SER SER ASN LEU PHE ASP LYS VAL SER PRO          
        VAL VAL SER HIS LYS VAL ASP LEU ALA VAL LEU ALA ALA          
        VAL GLU ILE ARG GLY VAL SER SER PRO ASP HIS VAL PHE          
        LEU PRO ILE PRO ASN TRP GLU HIS LYS GLU ASN PRO GLU          
        THR GLU GLU ASP VAL GLY PRO VAL VAL GLN HIS ILE TYR          
        GLU LEU ARG ASN ASN GLY PRO SER SER PHE SER LYS ALA          
        MET LEU HIS LEU GLN TRP PRO TYR LYS TYR ASN ASN ASN          
        THR LEU LEU TYR ILE LEU HIS TYR ASP ILE ASP GLY PRO          
        MET ASN CYS THR SER ASP MET GLU ILE ASN PRO LEU ARG          
        ILE LYS ILE SER SER LEU GLN THR THR GLU LYS ASN ASP          
        THR VAL ALA GLY GLN GLY GLU ARG ASP HIS LEU ILE THR          
        LYS ARG ASP LEU ALA LEU SER GLU GLY ASP ILE HIS THR          
        LEU GLY CYS GLY VAL ALA GLN CYS LEU LYS ILE VAL CYS          
        GLN VAL GLY ARG LEU ASP ARG GLY LYS SER ALA ILE LEU          
        TYR VAL LYS SER LEU LEU TRP THR GLU THR PHE MET ASN          
        LYS GLU ASN GLN ASN HIS SER TYR SER LEU LYS SER SER          
        ALA SER PHE ASN VAL ILE GLU PHE PRO TYR LYS ASN LEU          
        PRO ILE GLU ASP ILE THR ASN SER THR LEU VAL THR THR          
        ASN VAL THR TRP GLY ILE GLN PRO ALA PRO MET PRO VAL          
        PRO VAL TRP VAL ILE'''
        #
        self.authSequenceB='''
        GLY PRO ASN ILE CYS THR THR ARG GLY VAL SER SER CYS          
        GLN GLN CYS LEU ALA VAL SER PRO MET CYS ALA TRP CYS          
        SER ASP GLU ALA LEU PRO LEU GLY SER PRO ARG CYS ASP          
        LEU LYS GLU ASN LEU LEU LYS ASP ASN CYS ALA PRO GLU          
        SER ILE GLU PHE PRO VAL SER GLU ALA ARG VAL LEU GLU          
        ASP ARG PRO LEU SER ASP LYS GLY SER GLY ASP SER SER          
        GLN VAL THR GLN VAL SER PRO GLN ARG ILE ALA LEU ARG          
        LEU ARG PRO ASP ASP SER LYS ASN PHE SER ILE GLN VAL          
        ARG GLN VAL GLU ASP TYR PRO VAL ASP ILE TYR TYR LEU          
        MET ASP LEU SER TYR SER MET LYS ASP ASP LEU TRP SER          
        ILE GLN ASN LEU GLY THR LYS LEU ALA THR GLN MET ARG          
        LYS LEU THR SER ASN LEU ARG ILE GLY PHE GLY ALA PHE          
        VAL ASP LYS PRO VAL SER PRO TYR MET TYR ILE SER PRO          
        PRO GLU ALA LEU GLU ASN PRO CYS TYR ASP MET LYS THR          
        THR CYS LEU PRO MET PHE GLY TYR LYS HIS VAL LEU THR          
        LEU THR ASP GLN VAL THR ARG PHE ASN GLU GLU VAL LYS          
        LYS GLN SER VAL SER ARG ASN ARG ASP ALA PRO GLU GLY          
        GLY PHE ASP ALA ILE MET GLN ALA THR VAL CYS ASP GLU          
        LYS ILE GLY TRP ARG ASN ASP ALA SER HIS LEU LEU VAL          
        PHE THR THR ASP ALA LYS THR HIS ILE ALA LEU ASP GLY          
        ARG LEU ALA GLY ILE VAL GLN PRO ASN ASP GLY GLN CYS          
        HIS VAL GLY SER ASP ASN HIS TYR SER ALA SER THR THR          
        MET ASP TYR PRO SER LEU GLY LEU MET THR GLU LYS LEU          
        SER GLN LYS ASN ILE ASN LEU ILE PHE ALA VAL THR GLU          
        ASN VAL VAL ASN LEU TYR GLN ASN TYR SER GLU LEU ILE          
        PRO GLY THR THR VAL GLY VAL LEU SER MET ASP SER SER          
        ASN VAL LEU GLN LEU ILE VAL ASP ALA TYR GLY LYS ILE          
        ARG SER LYS VAL GLU LEU GLU VAL ARG ASP LEU PRO GLU          
        GLU LEU SER LEU SER PHE ASN ALA THR CYS LEU ASN ASN          
        GLU VAL ILE PRO GLY LEU LYS SER CYS MET GLY LEU LYS          
        ILE GLY ASP THR VAL SER PHE SER ILE GLU ALA LYS VAL          
        ARG GLY CYS PRO GLN GLU LYS GLU LYS SER PHE THR ILE          
        LYS PRO VAL GLY PHE LYS ASP SER LEU ILE VAL GLN VAL          
        THR PHE ASP CYS ASP CYS ALA CYS GLN ALA GLN ALA GLU          
        PRO ASN SER HIS ARG CYS ASN ASN GLY ASN GLY THR PHE          
        GLU CYS GLY VAL CYS ARG CYS GLY PRO GLY TRP LEU GLY          
        SER GLN CYS GLU CYS SER GLU GLU ASP TYR ARG PRO SER          
        GLN GLN ASP GLU CYS SER PRO ARG GLU GLY GLN PRO VAL          
        CYS SER GLN ARG GLY GLU CYS LEU CYS GLY GLN CYS VAL          
        CYS HIS SER SER ASP PHE GLY LYS ILE THR GLY LYS TYR          
        CYS GLU CYS ASP ASP PHE SER CYS VAL ARG TYR LYS GLY          
        GLU MET CYS SER GLY HIS GLY GLN CYS SER CYS GLY ASP          
        CYS LEU CYS ASP SER ASP TRP THR GLY TYR TYR CYS ASN          
        CYS THR THR ARG THR ASP THR CYS MET SER SER ASN GLY          
        LEU LEU CYS SER GLY ARG GLY LYS CYS GLU CYS GLY SER          
        CYS VAL CYS ILE GLN PRO GLY SER TYR GLY ASP THR CYS          
        GLU LYS CYS PRO THR CYS PRO ASP ALA CYS THR PHE LYS          
        LYS GLU CYS VAL GLU CYS LYS LYS PHE ASP ARG GLY ALA          
        LEU HIS ASP GLU ASN THR CYS ASN ARG TYR CYS ARG ASP          
        GLU ILE GLU SER VAL LYS GLU LEU LYS ASP THR GLY LYS          
        ASP ALA VAL ASN CYS THR TYR LYS ASN GLU ASP ASP CYS          
        VAL VAL ARG PHE GLN TYR TYR GLU ASP SER SER GLY LYS          
        SER ILE LEU TYR VAL VAL GLU GLU PRO GLU CYS PRO LYS          
        GLY PRO ASP ILE LEU VAL'''                                      

        self.xyzSequenceA=[
            ('PHE','1'), ('ASN','2'), ('LEU','3'), ('ASP','4'), ('VAL','5'),
            ('ASP','6'), ('SER','7'), ('PRO','8'), ('ALA','9'), ('GLU','10'), ('TYR','11'),
            ('SER','12'), ('GLY','13'), ('PRO','14'), ('GLU','15'), ('GLY','16'),
            ('SER','17'), ('TYR','18'), ('PHE','19'), ('GLY','20'), ('PHE','21'),
            ('ALA','22'), ('VAL','23'), ('ASP','24'), ('PHE','25'), ('PHE','26'),
            ('VAL','27'), ('PRO','28'), ('SER','29'), ('ALA','30'), ('SER','31'),
            ('SER','32'), ('ARG','33'), ('MET','34'), ('PHE','35'), ('LEU','36'),
            ('LEU','37'), ('VAL','38'), ('GLY','39'), ('ALA','40'), ('PRO','41'),
            ('LYS','42'), ('ALA','43'), ('ASN','44'), ('THR','45'), ('THR','46'),
            ('GLN','47'), ('PRO','48'), ('GLY','49'), ('ILE','50'), ('VAL','51'),
            ('GLU','52'), ('GLY','53'), ('GLY','54'), ('GLN','55'), ('VAL','56'),
            ('LEU','57'), ('LYS','58'), ('CYS','59'), ('ASP','60'), ('TRP','61'),
            ('SER','62'), ('SER','63'), ('THR','64'), ('ARG','65'), ('ARG','66'),
            ('CYS','67'), ('GLN','68'), ('PRO','69'), ('ILE','70'), ('GLU','71'),
            ('PHE','72'), ('ASP','73'), ('ALA','74'), ('THR','75'), ('GLY','76'),
            ('ASN','77'), ('ARG','78'), ('ASP','79'), ('TYR','80'), ('ALA','81'),
            ('LYS','82'), ('ASP','83'), ('ASP','84'), ('PRO','85'), ('LEU','86'),
            ('GLU','87'), ('PHE','88'), ('LYS','89'), ('SER','90'), ('HIS','91'),
            ('GLN','92'), ('TRP','93'), ('PHE','94'), ('GLY','95'), ('ALA','96'),
            ('SER','97'), ('VAL','98'), ('ARG','99'), ('SER','100'), ('LYS','101'),
            ('GLN','102'), ('ASP','103'), ('LYS','104'), ('ILE','105'), ('LEU','106'),
            ('ALA','107'), ('CYS','108'), ('ALA','109'), ('PRO','110'), ('LEU','111'),
            ('TYR','112'), ('HIS','113'), ('TRP','114'), ('ARG','115'), ('THR','116'),
            ('GLU','117'), ('MET','118'), ('LYS','119'), ('GLN','120'), ('GLU','121'),
            ('ARG','122'), ('GLU','123'), ('PRO','124'), ('VAL','125'), ('GLY','126'),
            ('THR','127'), ('CYS','128'), ('PHE','129'), ('LEU','130'), ('GLN','131'),
            ('ASP','132'), ('GLY','133'), ('THR','134'), ('LYS','135'), ('THR','136'),
            ('VAL','137'), ('GLU','138'), ('TYR','139'), ('ALA','140'), ('PRO','141'),
            ('CYS','142'), ('ARG','143'), ('SER','144'), ('GLN','145'), ('ASP','146'),
            ('ILE','147'), ('ASP','148'), ('ALA','149'), ('ASP','150'), ('GLY','151'),
            ('GLN','152'), ('GLY','153'), ('PHE','154'), ('CYS','155'), ('GLN','156'),
            ('GLY','157'), ('GLY','158'), ('PHE','159'), ('SER','160'), ('ILE','161'),
            ('ASP','162'), ('PHE','163'), ('THR','164'), ('LYS','165'), ('ALA','166'),
            ('ASP','167'), ('ARG','168'), ('VAL','169'), ('LEU','170'), ('LEU','171'),
            ('GLY','172'), ('GLY','173'), ('PRO','174'), ('GLY','175'), ('SER','176'),
            ('PHE','177'), ('TYR','178'), ('TRP','179'), ('GLN','180'), ('GLY','181'),
            ('GLN','182'), ('LEU','183'), ('ILE','184'), ('SER','185'), ('ASP','186'),
            ('GLN','187'), ('VAL','188'), ('ALA','189'), ('GLU','190'), ('ILE','191'),
            ('VAL','192'), ('SER','193'), ('LYS','194'), ('TYR','195'), ('ASP','196'),
            ('PRO','197'), ('ASN','198'), ('VAL','199'), ('TYR','200'), ('SER','201'),
            ('ILE','202'), ('LYS','203'), ('TYR','204'), ('ASN','205'), ('ASN','206'),
            ('GLN','207'), ('LEU','208'), ('ALA','209'), ('THR','210'), ('ARG','211'),
            ('THR','212'), ('ALA','213'), ('GLN','214'), ('ALA','215'), ('ILE','216'),
            ('PHE','217'), ('ASP','218'), ('ASP','219'), ('SER','220'), ('TYR','221'),
            ('LEU','222'), ('GLY','223'), ('TYR','224'), ('SER','225'), ('VAL','226'),
            ('ALA','227'), ('VAL','228'), ('GLY','229'), ('ASP','230'), ('PHE','231'),
            ('ASN','232'), ('GLY','233'), ('ASP','234'), ('GLY','235'), ('ILE','236'),
            ('ASP','237'), ('ASP','238'), ('PHE','239'), ('VAL','240'), ('SER','241'),
            ('GLY','242'), ('VAL','243'), ('PRO','244'), ('ARG','245'), ('ALA','246'),
            ('ALA','247'), ('ARG','248'), ('THR','249'), ('LEU','250'), ('GLY','251'),
            ('MET','252'), ('VAL','253'), ('TYR','254'), ('ILE','255'), ('TYR','256'),
            ('ASP','257'), ('GLY','258'), ('LYS','259'), ('ASN','260'), ('MET','261'),
            ('SER','262'), ('SER','263'), ('LEU','264'), ('TYR','265'), ('ASN','266'),
            ('PHE','267'), ('THR','268'), ('GLY','269'), ('GLU','270'), ('GLN','271'),
            ('MET','272'), ('ALA','273'), ('ALA','274'), ('TYR','275'), ('PHE','276'),
            ('GLY','277'), ('PHE','278'), ('SER','279'), ('VAL','280'), ('ALA','281'),
            ('ALA','282'), ('THR','283'), ('ASP','284'), ('ILE','285'), ('ASN','286'),
            ('GLY','287'), ('ASP','288'), ('ASP','289'), ('TYR','290'), ('ALA','291'),
            ('ASP','292'), ('VAL','293'), ('PHE','294'), ('ILE','295'), ('GLY','296'),
            ('ALA','297'), ('PRO','298'), ('LEU','299'), ('PHE','300'), ('MET','301'),
            ('ASP','302'), ('ARG','303'), ('GLY','304'), ('SER','305'), ('ASP','306'),
            ('GLY','307'), ('LYS','308'), ('LEU','309'), ('GLN','310'), ('GLU','311'),
            ('VAL','312'), ('GLY','313'), ('GLN','314'), ('VAL','315'), ('SER','316'),
            ('VAL','317'), ('SER','318'), ('LEU','319'), ('GLN','320'), ('ARG','321'),
            ('ALA','322'), ('SER','323'), ('GLY','324'), ('ASP','325'), ('PHE','326'),
            ('GLN','327'), ('THR','328'), ('THR','329'), ('LYS','330'), ('LEU','331'),
            ('ASN','332'), ('GLY','333'), ('PHE','334'), ('GLU','335'), ('VAL','336'),
            ('PHE','337'), ('ALA','338'), ('ARG','339'), ('PHE','340'), ('GLY','341'),
            ('SER','342'), ('ALA','343'), ('ILE','344'), ('ALA','345'), ('PRO','346'),
            ('LEU','347'), ('GLY','348'), ('ASP','349'), ('LEU','350'), ('ASP','351'),
            ('GLN','352'), ('ASP','353'), ('GLY','354'), ('PHE','355'), ('ASN','356'),
            ('ASP','357'), ('ILE','358'), ('ALA','359'), ('ILE','360'), ('ALA','361'),
            ('ALA','362'), ('PRO','363'), ('TYR','364'), ('GLY','365'), ('GLY','366'),
            ('GLU','367'), ('ASP','368'), ('LYS','369'), ('LYS','370'), ('GLY','371'),
            ('ILE','372'), ('VAL','373'), ('TYR','374'), ('ILE','375'), ('PHE','376'),
            ('ASN','377'), ('GLY','378'), ('ARG','379'), ('SER','380'), ('THR','381'),
            ('GLY','382'), ('LEU','383'), ('ASN','384'), ('ALA','385'), ('VAL','386'),
            ('PRO','387'), ('SER','388'), ('GLN','389'), ('ILE','390'), ('LEU','391'),
            ('GLU','392'), ('GLY','393'), ('GLN','394'), ('TRP','395'), ('ALA','396'),
            ('ALA','397'), ('ARG','398'), ('SER','399'), ('MET','400'), ('PRO','401'),
            ('PRO','402'), ('SER','403'), ('PHE','404'), ('GLY','405'), ('TYR','406'),
            ('SER','407'), ('MET','408'), ('LYS','409'), ('GLY','410'), ('ALA','411'),
            ('THR','412'), ('ASP','413'), ('ILE','414'), ('ASP','415'), ('LYS','416'),
            ('ASN','417'), ('GLY','418'), ('TYR','419'), ('PRO','420'), ('ASP','421'),
            ('LEU','422'), ('ILE','423'), ('VAL','424'), ('GLY','425'), ('ALA','426'),
            ('PHE','427'), ('GLY','428'), ('VAL','429'), ('ASP','430'), ('ARG','431'),
            ('ALA','432'), ('ILE','433'), ('LEU','434'), ('TYR','435'), ('ARG','436'),
            ('ALA','437'), ('ARG','438'), ('PRO','439'), ('VAL','440'), ('ILE','441'),
            ('THR','442'), ('VAL','443'), ('ASN','444'), ('ALA','445'), ('GLY','446'),
            ('LEU','447'), ('GLU','448'), ('VAL','449'), ('TYR','450'), ('PRO','451'),
            ('SER','452'), ('ILE','453'), ('LEU','454'), ('ASN','455'), ('GLN','456'),
            ('ASP','457'), ('ASN','458'), ('LYS','459'), ('THR','460'), ('CYS','461'),
            ('SER','462'), ('LEU','463'), ('PRO','464'), ('GLY','465'), ('THR','466'),
            ('ALA','467'), ('LEU','468'), ('LYS','469'), ('VAL','470'), ('SER','471'),
            ('CYS','472'), ('PHE','473'), ('ASN','474'), ('VAL','475'), ('ARG','476'),
            ('PHE','477'), ('CYS','478'), ('LEU','479'), ('LYS','480'), ('ALA','481'),
            ('ASP','482'), ('GLY','483'), ('LYS','484'), ('GLY','485'), ('VAL','486'),
            ('LEU','487'), ('PRO','488'), ('ARG','489'), ('LYS','490'), ('LEU','491'),
            ('ASN','492'), ('PHE','493'), ('GLN','494'), ('VAL','495'), ('GLU','496'),
            ('LEU','497'), ('LEU','498'), ('LEU','499'), ('ASP','500'), ('LYS','501'),
            ('LEU','502'), ('LYS','503'), ('GLN','504'), ('LYS','505'), ('GLY','506'),
            ('ALA','507'), ('ILE','508'), ('ARG','509'), ('ARG','510'), ('ALA','511'),
            ('LEU','512'), ('PHE','513'), ('LEU','514'), ('TYR','515'), ('SER','516'),
            ('ARG','517'), ('SER','518'), ('PRO','519'), ('SER','520'), ('HIS','521'),
            ('SER','522'), ('LYS','523'), ('ASN','524'), ('MET','525'), ('THR','526'),
            ('ILE','527'), ('SER','528'), ('ARG','529'), ('GLY','530'), ('GLY','531'),
            ('LEU','532'), ('MET','533'), ('GLN','534'), ('CYS','535'), ('GLU','536'),
            ('GLU','537'), ('LEU','538'), ('ILE','539'), ('ALA','540'), ('TYR','541'),
            ('LEU','542'), ('ARG','543'), ('ASP','544'), ('GLU','545'), ('SER','546'),
            ('GLU','547'), ('PHE','548'), ('ARG','549'), ('ASP','550'), ('LYS','551'),
            ('LEU','552'), ('THR','553'), ('PRO','554'), ('ILE','555'), ('THR','556'),
            ('ILE','557'), ('PHE','558'), ('MET','559'), ('GLU','560'), ('TYR','561'),
            ('ARG','562'), ('LEU','563'), ('ASP','564'), ('TYR','565'), ('ARG','566'),
            ('THR','567'), ('ALA','568'), ('ALA','569'), ('ASP','570'), ('THR','571'),
            ('THR','572'), ('GLY','573'), ('LEU','574'), ('GLN','575'), ('PRO','576'),
            ('ILE','577'), ('LEU','578'), ('ASN','579'), ('GLN','580'), ('PHE','581'),
            ('THR','582'), ('PRO','583'), ('ALA','584'), ('ASN','585'), ('ILE','586'),
            ('SER','587'), ('ARG','588'), ('GLN','589'), ('ALA','590'), ('HIS','591'),
            ('ILE','592'), ('LEU','593'), ('LEU','594'), ('ASP','595'), ('CYS','596'),
            ('GLY','597'), ('GLU','598'), ('ASP','599'), ('ASN','600'), ('VAL','601'),
            ('CYS','602'), ('LYS','603'), ('PRO','604'), ('LYS','605'), ('LEU','606'),
            ('GLU','607'), ('VAL','608'), ('SER','609'), ('VAL','610'), ('ASP','611'),
            ('SER','612'), ('ASP','613'), ('GLN','614'), ('LYS','615'), ('LYS','616'),
            ('ILE','617'), ('TYR','618'), ('ILE','619'), ('GLY','620'), ('ASP','621'),
            ('ASP','622'), ('ASN','623'), ('PRO','624'), ('LEU','625'), ('THR','626'),
            ('LEU','627'), ('ILE','628'), ('VAL','629'), ('LYS','630'), ('ALA','631'),
            ('GLN','632'), ('ASN','633'), ('GLN','634'), ('GLY','635'), ('GLU','636'),
            ('GLY','637'), ('ALA','638'), ('TYR','639'), ('GLU','640'), ('ALA','641'),
            ('GLU','642'), ('LEU','643'), ('ILE','644'), ('VAL','645'), ('SER','646'),
            ('ILE','647'), ('PRO','648'), ('LEU','649'), ('GLN','650'), ('ALA','651'),
            ('ASP','652'), ('PHE','653'), ('ILE','654'), ('GLY','655'), ('VAL','656'),
            ('VAL','657'), ('ARG','658'), ('ASN','659'), ('ASN','660'), ('GLU','661'),
            ('ALA','662'), ('LEU','663'), ('ALA','664'), ('ARG','665'), ('LEU','666'),
            ('SER','667'), ('CYS','668'), ('ALA','669'), ('PHE','670'), ('LYS','671'),
            ('THR','672'), ('GLU','673'), ('ASN','674'), ('GLN','675'), ('THR','676'),
            ('ARG','677'), ('GLN','678'), ('VAL','679'), ('VAL','680'), ('CYS','681'),
            ('ASP','682'), ('LEU','683'), ('GLY','684'), ('ASN','685'), ('PRO','686'),
            ('MET','687'), ('LYS','688'), ('ALA','689'), ('GLY','690'), ('THR','691'),
            ('GLN','692'), ('LEU','693'), ('LEU','694'), ('ALA','695'), ('GLY','696'),
            ('LEU','697'), ('ARG','698'), ('PHE','699'), ('SER','700'), ('VAL','701'),
            ('HIS','702'), ('GLN','703'), ('GLN','704'), ('SER','705'), ('GLU','706'),
            ('MET','707'), ('ASP','708'), ('THR','709'), ('SER','710'), ('VAL','711'),
            ('LYS','712'), ('PHE','713'), ('ASP','714'), ('LEU','715'), ('GLN','716'),
            ('ILE','717'), ('GLN','718'), ('SER','719'), ('SER','720'), ('ASN','721'),
            ('LEU','722'), ('PHE','723'), ('ASP','724'), ('LYS','725'), ('VAL','726'),
            ('SER','727'), ('PRO','728'), ('VAL','729'), ('VAL','730'), ('SER','731'),
            ('HIS','732'), ('LYS','733'), ('VAL','734'), ('ASP','735'), ('LEU','736'),
            ('ALA','737'), ('VAL','738'), ('LEU','739'), ('ALA','740'), ('ALA','741'),
            ('VAL','742'), ('GLU','743'), ('ILE','744'), ('ARG','745'), ('GLY','746'),
            ('VAL','747'), ('SER','748'), ('SER','749'), ('PRO','750'), ('ASP','751'),
            ('HIS','752'), ('VAL','753'), ('PHE','754'), ('LEU','755'), ('PRO','756'),
            ('ILE','757'), ('PRO','758'), ('ASN','759'), ('TRP','760'), ('GLU','761'),
            ('HIS','762'), ('LYS','763'), ('GLU','764'), ('ASN','765'), ('PRO','766'),
            ('GLU','767'), ('THR','768'), ('GLU','769'), ('GLU','770'), ('ASP','771'),
            ('VAL','772'), ('GLY','773'), ('PRO','774'), ('VAL','775'), ('VAL','776'),
            ('GLN','777'), ('HIS','778'), ('ILE','779'), ('TYR','780'), ('GLU','781'),
            ('LEU','782'), ('ARG','783'), ('ASN','784'), ('ASN','785'), ('GLY','786'),
            ('PRO','787'), ('SER','788'), ('SER','789'), ('PHE','790'), ('SER','791'),
            ('LYS','792'), ('ALA','793'), ('MET','794'), ('LEU','795'), ('HIS','796'),
            ('LEU','797'), ('GLN','798'), ('TRP','799'), ('PRO','800'), ('TYR','801'),
            ('LYS','802'), ('TYR','803'), ('ASN','804'), ('ASN','805'), ('ASN','806'),
            ('THR','807'), ('LEU','808'), ('LEU','809'), ('TYR','810'), ('ILE','811'),
            ('LEU','812'), ('HIS','813'), ('TYR','814'), ('ASP','815'), ('ILE','816'),
            ('ASP','817'), ('GLY','818'), ('PRO','819'), ('MET','820'), ('ASN','821'),
            ('CYS','822'), ('THR','823'), ('SER','824'), ('ASP','825'), ('MET','826'),
            ('GLU','827'), ('ILE','828'), ('ASN','829'), ('PRO','830'), ('LEU','831'),
            ('ARG','832'), ('ILE','833'), ('LYS','834'), ('ILE','835'), ('SER','836'),
            ('SER','837'), ('LEU','838'), ('ASP','868'), ('ILE','869'), ('HIS','870'),
            ('THR','871'), ('LEU','872'), ('GLY','873'), ('CYS','874'), ('GLY','875'),
            ('VAL','876'), ('ALA','877'), ('GLN','878'), ('CYS','879'), ('LEU','880'),
            ('LYS','881'), ('ILE','882'), ('VAL','883'), ('CYS','884'), ('GLN','885'),
            ('VAL','886'), ('GLY','887'), ('ARG','888'), ('LEU','889'), ('ASP','890'),
            ('ARG','891'), ('GLY','892'), ('LYS','893'), ('SER','894'), ('ALA','895'),
            ('ILE','896'), ('LEU','897'), ('TYR','898'), ('VAL','899'), ('LYS','900'),
            ('SER','901'), ('LEU','902'), ('LEU','903'), ('TRP','904'), ('THR','905'),
            ('GLU','906'), ('THR','907'), ('PHE','908'), ('MET','909'), ('ASN','910'),
            ('LYS','911'), ('GLU','912'), ('ASN','913'), ('GLN','914'), ('ASN','915'),
            ('HIS','916'), ('SER','917'), ('TYR','918'), ('SER','919'), ('LEU','920'),
            ('LYS','921'), ('SER','922'), ('SER','923'), ('ALA','924'), ('SER','925'),
            ('PHE','926'), ('ASN','927'), ('VAL','928'), ('ILE','929'), ('GLU','930'),
            ('PHE','931'), ('PRO','932'), ('TYR','933'), ('LYS','934'), ('ASN','935'),
            ('LEU','936'), ('PRO','937'), ('ILE','938'), ('GLU','939'), ('ASP','940'),
            ('ILE','941'), ('THR','942'), ('ASN','943'), ('SER','944'), ('THR','945'),
            ('LEU','946'), ('VAL','947'), ('THR','948'), ('THR','949'), ('ASN','950'),
            ('VAL','951'), ('THR','952'), ('TRP','953'), ('GLY','954'), ('ILE','955'),
            ('GLN','956'), ('PRO','957'), ('ALA','958'), ('PRO','959'), ('MET','960'),
            ('PRO','961'), ('VAL','962'), ('PRO','963'), ('VAL','964'), ('TRP','965'),
            ('VAL','966'), ('ILE','967')
            ] 
        self.xyzSequenceAMod=[
            ('PHE','1'), ('ASN','2'), ('LEU','3'), ('ASP','4'), ('VAL','5'),
            ('ASP','6'), ('SER','7'), ('PRO','8'), ('ALA','9'), ('GLU','10'), ('TYR','11'),
            ('SER','12'), ('ALA','13'), ('PRO','14'), ('GLU','15'), ('ALA','16'),
            ('SER','17'), ('TYR','18'), ('PHE','19'), ('ALA','20'), ('PHE','21'),
            ('ALA','22'), ('VAL','23'), ('ASP','24'), ('PHE','25'), ('PHE','26'),
            ('VAL','27'), ('PRO','28'), ('SER','29'), ('ALA','30'), ('SER','31'),
            ('SER','32'), ('ARG','33'), ('MET','34'), ('PHE','35'), ('LEU','36'),
            ('LEU','37'), ('VAL','38'), ('ALA','39'), ('ALA','40'), ('PRO','41'),
            ('LYS','42'), ('ALA','43'), ('ASN','44'), ('THR','45'), ('THR','46'),
            ('GLN','47'), ('PRO','48'), ('ALA','49'), ('ILE','50'), ('VAL','51'),
            ('GLU','52'), ('ALA','53'), ('ALA','54'), ('GLN','55'), ('VAL','56'),
            ('LEU','57'), ('LYS','58'), ('CYS','59'), ('ASP','60'), ('TRP','61'),
            ('SER','62'), ('SER','63'), ('THR','64'), ('ARG','65'), ('ARG','66'),
            ('CYS','67'), ('GLN','68'), ('PRO','69'), ('ILE','70'), ('GLU','71'),
            ('PHE','72'), ('ASP','73'), ('ALA','74'), ('THR','75'), ('ALA','76'),
            ('ASN','77'), ('ARG','78'), ('ASP','79'), ('TYR','80'), ('ALA','81'),
            ('LYS','82'), ('ASP','83'), ('ASP','84'), ('PRO','85'), ('LEU','86'),
            ('GLU','87'), ('PHE','88'), ('LYS','89'), ('SER','90'), ('HIS','91'),
            ('GLN','92'), ('TRP','93'), ('PHE','94'), ('ALA','95'), ('ALA','96'),
            ('SER','97'), ('VAL','98'), ('ARG','99'), ('SER','100'), ('LYS','101'),
            ('GLN','102'), ('ASP','103'), ('LYS','104'), ('ILE','105'), ('LEU','106'),
            ('ALA','107'), ('CYS','108'), ('ALA','109'), ('PRO','110'), ('LEU','111'),
            ('TYR','112'), ('HIS','113'), ('TRP','114'), ('ARG','115'), ('THR','116'),
            ('GLU','117'), ('MET','118'), ('LYS','119'), ('GLN','120'), ('GLU','121'),
            ('ARG','122'), ('GLU','123'), ('PRO','124'), ('VAL','125'), ('ALA','126'),
            ('THR','127'), ('CYS','128'), ('PHE','129'), ('LEU','130'), ('GLN','131'),
            ('ASP','132'), ('ALA','133'), ('THR','134'), ('LYS','135'), ('THR','136'),
            ('VAL','137'), ('GLU','138'), ('TYR','139'), ('ALA','140'), ('PRO','141'),
            ('CYS','142'), ('ARG','143'), ('SER','144'), ('GLN','145'), ('ASP','146'),
            ('ILE','147'), ('ASP','148'), ('ALA','149'), ('ASP','150'), ('ALA','151'),
            ('GLN','152'), ('ALA','153'), ('PHE','154'), ('CYS','155'), ('GLN','156'),
            ('ALA','157'), ('ALA','158'), ('PHE','159'), ('SER','160'), ('ILE','161'),
            ('ASP','162'), ('PHE','163'), ('THR','164'), ('LYS','165'), ('ALA','166'),
            ('ASP','167'), ('ARG','168'), ('VAL','169'), ('LEU','170'), ('LEU','171'),
            ('ALA','172'), ('ALA','173'), ('PRO','174'), ('ALA','175'), ('SER','176'),
            ('PHE','177'), ('TYR','178'), ('TRP','179'), ('GLN','180'), ('ALA','181'),
            ('GLN','182'), ('LEU','183'), ('ILE','184'), ('SER','185'), ('ASP','186'),
            ('GLN','187'), ('VAL','188'), ('ALA','189'), ('GLU','190'), ('ILE','191'),
            ('VAL','192'), ('SER','193'), ('LYS','194'), ('TYR','195'), ('ASP','196'),
            ('PRO','197'), ('ASN','198'), ('VAL','199'), ('TYR','200'), ('SER','201'),
            ('ILE','202'), ('LYS','203'), ('TYR','204'), ('ASN','205'), ('ASN','206'),
            ('GLN','207'), ('LEU','208'), ('ALA','209'), ('THR','210'), ('ARG','211'),
            ('THR','212'), ('ALA','213'), ('GLN','214'), ('ALA','215'), ('ILE','216'),
            ('PHE','217'), ('ASN','218'), ('ASN','219'), ('SER','220'), ('TYR','221'),
            ('LEU','222'), ('ALA','223'), ('TYR','224'), ('SER','225'), ('VAL','226'),
            ('ALA','227'), ('VAL','228'), ('ALA','229'), ('ASN','230'), ('PHE','231'),
            ('ASN','232'), ('ALA','233'), ('ASN','234'), ('ALA','235'), ('ILE','236'),
            ('ASN','237'), ('ASN','238'), ('PHE','239'), ('VAL','240'), ('SER','241'),
            ('ALA','242'), ('VAL','243'), ('PRO','244'), ('ARG','245'), ('ALA','246'),
            ('ALA','247'), ('ARG','248'), ('THR','249'), ('LEU','250'), ('ALA','251'),
            ('MET','252'), ('VAL','253'), ('TYR','254'), ('ILE','255'), ('TYR','256'),
            ('ASN','257'), ('ALA','258'), ('LYS','259'), ('ASN','260'), ('MET','261'),
            ('SER','262'), ('SER','263'), ('LEU','264'), ('TYR','265'), ('ASN','266'),
            ('PHE','267'), ('THR','268'), ('ALA','269'), ('GLU','270'), ('GLN','271'),
            ('MET','272'), ('ALA','273'), ('ALA','274'), ('TYR','275'), ('PHE','276'),
            ('ALA','277'), ('PHE','278'), ('SER','279'), ('VAL','280'), ('ALA','281'),
            ('ALA','282'), ('THR','283'), ('ASN','284'), ('ILE','285'), ('ASN','286'),
            ('ALA','287'), ('ASN','288'), ('ASN','289'), ('TYR','290'), ('ALA','291'),
            ('ASP','292'), ('VAL','293'), ('PHE','294'), ('ILE','295'), ('ALA','296'),
            ('ALA','297'), ('PRO','298'), ('LEU','299'), ('PHE','300'), ('MET','301'),
            ('ASP','302'), ('ARG','303'), ('ALA','304'), ('SER','305'), ('ASP','306'),
            ('ALA','307'), ('LYS','308'), ('LEU','309'), ('GLN','310'), ('GLN','311'),
            ('VAL','312'), ('ALA','313'), ('GLN','314'), ('VAL','315'), ('SER','316'),
            ('VAL','317'), ('SER','318'), ('LEU','319'), ('GLN','320'), ('ARG','321'),
            ('ALA','322'), ('SER','323'), ('ALA','324'), ('ASP','325'), ('PHE','326'),
            ('GLN','327'), ('THR','328'), ('THR','329'), ('LYS','330'), ('LEU','331'),
            ('ASN','332'), ('GLY','333'), ('PHE','334'), ('GLN','335'), ('VAL','336'),
            ('PHE','337'), ('ALA','338'), ('ARG','339'), ('PHE','340'), ('GLY','341'),
            ('SER','342'), ('ALA','343'), ('ILE','344'), ('ALA','345'), ('PRO','346'),
            ('LEU','347'), ('GLY','348'), ('ASP','349'), ('LEU','350'), ('ASP','351'),
            ('GLN','352'), ('ASP','353'), ('GLY','354'), ('PHE','355'), ('ASN','356'),
            ('ASP','357'), ('ILE','358'), ('ALA','359'), ('ILE','360'), ('ALA','361'),
            ('ALA','362'), ('PRO','363'), ('TYR','364'), ('GLY','365'), ('GLY','366'),
            ('GLN','367'), ('ASP','368'), ('LYS','369'), ('LYS','370'), ('GLY','371'),
            ('ILE','372'), ('VAL','373'), ('TYR','374'), ('ILE','375'), ('PHE','376'),
            ('ASN','377'), ('GLY','378'), ('ARG','379'), ('SER','380'), ('THR','381'),
            ('GLY','382'), ('LEU','383'), ('ASN','384'), ('ALA','385'), ('VAL','386'),
            ('PRO','387'), ('SER','388'), ('GLN','389'), ('ILE','390'), ('LEU','391'),
            ('GLN','392'), ('GLY','393'), ('GLN','394'), ('TRP','395'), ('ALA','396'),
            ('ALA','397'), ('ARG','398'), ('SER','399'), ('MET','400'), ('PRO','401'),
            ('PRO','402'), ('SER','403'), ('PHE','404'), ('GLY','405'), ('TYR','406'),
            ('SER','407'), ('MET','408'), ('LYS','409'), ('GLY','410'), ('ALA','411'),
            ('THR','412'), ('ASP','413'), ('ILE','414'), ('ASP','415'), ('LYS','416'),
            ('ASN','417'), ('GLY','418'), ('TYR','419'), ('PRO','420'), ('ASP','421'),
            ('LEU','422'), ('ILE','423'), ('VAL','424'), ('GLY','425'), ('ALA','426'),
            ('PHE','427'), ('GLY','428'), ('VAL','429'), ('ASP','430'), ('ARG','431'),
            ('ALA','432'), ('ILE','433'), ('LEU','434'), ('TYR','435'), ('ARG','436'),
            ('ALA','437'), ('ARG','438'), ('PRO','439'), ('VAL','440'), ('ILE','441'),
            ('THR','442'), ('VAL','443'), ('ASN','444'), ('ALA','445'), ('GLY','446'),
            ('LEU','447'), ('GLN','448'), ('VAL','449'), ('TYR','450'), ('PRO','451'),
            ('SER','452'), ('ILE','453'), ('LEU','454'), ('ASN','455'), ('GLN','456'),
            ('ASP','457'), ('ASN','458'), ('LYS','459'), ('THR','460'), ('CYS','461'),
            ('SER','462'), ('LEU','463'), ('PRO','464'), ('GLY','465'), ('THR','466'),
            ('ALA','467'), ('LEU','468'), ('LYS','469'), ('VAL','470'), ('SER','471'),
            ('CYS','472'), ('PHE','473'), ('ASN','474'), ('VAL','475'), ('ARG','476'),
            ('PHE','477'), ('CYS','478'), ('LEU','479'), ('LYS','480'), ('ALA','481'),
            ('ASP','482'), ('GLY','483'), ('LYS','484'), ('GLY','485'), ('VAL','486'),
            ('LEU','487'), ('PRO','488'), ('ARG','489'), ('LYS','490'), ('LEU','491'),
            ('ASN','492'), ('PHE','493'), ('GLN','494'), ('VAL','495'), ('GLN','496'),
            ('LEU','497'), ('LEU','498'), ('LEU','499'), ('ASP','500'), ('LYS','501'),
            ('LEU','502'), ('LYS','503'), ('GLN','504'), ('LYS','505'), ('GLY','506'),
            ('ALA','507'), ('ILE','508'), ('ARG','509'), ('ARG','510'), ('ALA','511'),
            ('LEU','512'), ('PHE','513'), ('LEU','514'), ('TYR','515'), ('SER','516'),
            ('ARG','517'), ('SER','518'), ('PRO','519'), ('SER','520'), ('HIS','521'),
            ('SER','522'), ('LYS','523'), ('ASN','524'), ('MET','525'), ('THR','526'),
            ('ILE','527'), ('SER','528'), ('ARG','529'), ('GLY','530'), ('GLY','531'),
            ('LEU','532'), ('MET','533'), ('GLN','534'), ('CYS','535'), ('GLN','536'),
            ('GLN','537'), ('LEU','538'), ('ILE','539'), ('ALA','540'), ('TYR','541'),
            ('LEU','542'), ('ARG','543'), ('ASP','544'), ('GLN','545'), ('SER','546'),
            ('GLN','547'), ('PHE','548'), ('ARG','549'), ('ASP','550'), ('LYS','551'),
            ('LEU','552'), ('THR','553'), ('PRO','554'), ('ILE','555'), ('THR','556'),
            ('ILE','557'), ('PHE','558'), ('MET','559'), ('GLN','560'), ('TYR','561'),
            ('ARG','562'), ('LEU','563'), ('ASP','564'), ('TYR','565'), ('ARG','566'),
            ('THR','567'), ('ALA','568'), ('ALA','569'), ('ASP','570'), ('THR','571'),
            ('THR','572'), ('GLY','573'), ('LEU','574'), ('GLN','575'), ('PRO','576'),
            ('ILE','577'), ('LEU','578'), ('ASN','579'), ('GLN','580'), ('PHE','581'),
            ('THR','582'), ('PRO','583'), ('ALA','584'), ('ASN','585'), ('ILE','586'),
            ('SER','587'), ('ARG','588'), ('GLN','589'), ('ALA','590'), ('HIS','591'),
            ('ILE','592'), ('LEU','593'), ('LEU','594'), ('ASP','595'), ('CYS','596'),
            ('GLY','597'), ('GLN','598'), ('ASP','599'), ('ASN','600'), ('VAL','601'),
            ('CYS','602'), ('LYS','603'), ('PRO','604'), ('LYS','605'), ('LEU','606'),
            ('GLU','607'), ('VAL','608'), ('SER','609'), ('VAL','610'), ('ASP','611'),
            ('SER','612'), ('ASP','613'), ('GLN','614'), ('LYS','615'), ('LYS','616'),
            ('ILE','617'), ('TYR','618'), ('ILE','619'), ('GLY','620'), ('ASP','621'),
            ('ASP','622'), ('ASN','623'), ('PRO','624'), ('LEU','625'), ('THR','626'),
            ('LEU','627'), ('ILE','628'), ('VAL','629'), ('LYS','630'), ('ALA','631'),
            ('GLN','632'), ('ASN','633'), ('GLN','634'), ('GLY','635'), ('GLU','636'),
            ('GLY','637'), ('ALA','638'), ('TYR','639'), ('GLU','640'), ('ALA','641'),
            ('GLU','642'), ('LEU','643'), ('ILE','644'), ('VAL','645'), ('SER','646'),
            ('ILE','647'), ('PRO','648'), ('LEU','649'), ('GLN','650'), ('ALA','651'),
            ('ASP','652'), ('PHE','653'), ('ILE','654'), ('GLY','655'), ('VAL','656'),
            ('VAL','657'), ('ARG','658'), ('ASN','659'), ('ASN','660'), ('GLU','661'),
            ('ALA','662'), ('LEU','663'), ('ALA','664'), ('ARG','665'), ('LEU','666'),
            ('SER','667'), ('CYS','668'), ('ALA','669'), ('PHE','670'), ('LYS','671'),
            ('THR','672'), ('GLU','673'), ('ASN','674'), ('GLN','675'), ('THR','676'),
            ('ARG','677'), ('GLN','678'), ('VAL','679'), ('VAL','680'), ('CYS','681'),
            ('ASP','682'), ('LEU','683'), ('GLY','684'), ('ASN','685'), ('PRO','686'),
            ('MET','687'), ('LYS','688'), ('ALA','689'), ('GLY','690'), ('THR','691'),
            ('GLN','692'), ('LEU','693'), ('LEU','694'), ('ALA','695'), ('GLY','696'),
            ('LEU','697'), ('ARG','698'), ('PHE','699'), ('SER','700'), ('VAL','701'),
            ('HIS','702'), ('GLN','703'), ('GLN','704'), ('SER','705'), ('GLU','706'),
            ('MET','707'), ('ASP','708'), ('THR','709'), ('SER','710'), ('VAL','711'),
            ('LYS','712'), ('PHE','713'), ('ASP','714'), ('LEU','715'), ('GLN','716'),
            ('ILE','717'), ('GLN','718'), ('SER','719'), ('SER','720'), ('ASN','721'),
            ('LEU','722'), ('PHE','723'), ('ASP','724'), ('LYS','725'), ('VAL','726'),
            ('SER','727'), ('PRO','728'), ('VAL','729'), ('VAL','730'), ('SER','731'),
            ('HIS','732'), ('LYS','733'), ('VAL','734'), ('ASN','735'), ('LEU','736'),
            ('ALA','737'), ('VAL','738'), ('LEU','739'), ('ALA','740'), ('ALA','741'),
            ('VAL','742'), ('GLU','743'), ('ILE','744'), ('ARG','745'), ('GLY','746'),
            ('VAL','747'), ('SER','748'), ('SER','749'), ('PRO','750'), ('ASN','751'),
            ('HIS','752'), ('VAL','753'), ('PHE','754'), ('LEU','755'), ('PRO','756'),
            ('ILE','757'), ('PRO','758'), ('ASN','759'), ('TRP','760'), ('GLU','761'),
            ('HIS','762'), ('LYS','763'), ('GLN','764'), ('ASN','765'), ('PRO','766'),
            ('GLN','767'), ('THR','768'), ('GLN','769'), ('GLN','770'), ('ASN','771'),
            ('VAL','772'), ('GLY','773'), ('PRO','774'), ('VAL','775'), ('VAL','776'),
            ('GLN','777'), ('HIS','778'), ('ILE','779'), ('TYR','780'), ('GLN','781'),
            ('LEU','782'), ('ARG','783'), ('ASN','784'), ('ASN','785'), ('GLY','786'),
            ('PRO','787'), ('SER','788'), ('SER','789'), ('PHE','790'), ('SER','791'),
            ('LYS','792'), ('ALA','793'), ('MET','794'), ('LEU','795'), ('HIS','796'),
            ('LEU','797'), ('GLN','798'), ('TRP','799'), ('PRO','800'), ('TYR','801'),
            ('LYS','802'), ('TYR','803'), ('ASN','804'), ('ASN','805'), ('ASN','806'),
            ('THR','807'), ('LEU','808'), ('LEU','809'), ('TYR','810'), ('ILE','811'),
            ('LEU','812'), ('HIS','813'), ('TYR','814'), ('ASN','815'), ('ILE','816'),
            ('ASN','817'), ('GLY','818'), ('PRO','819'), ('MET','820'), ('ASN','821'),
            ('CYS','822'), ('THR','823'), ('SER','824'), ('ASN','825'), ('MET','826'),
            ('GLN','827'), ('ILE','828'), ('ASN','829'), ('PRO','830'), ('LEU','831'),
            ('ARG','832'), ('ILE','833'), ('LYS','834'), ('ILE','835'), ('SER','836'),
            ('SER','837'), ('LEU','838'), ('ASN','868'), ('ILE','869'), ('HIS','870'),
            ('THR','871'), ('LEU','872'), ('GLY','873'), ('CYS','874'), ('GLY','875'),
            ('VAL','876'), ('ALA','877'), ('GLN','878'), ('CYS','879'), ('LEU','880'),
            ('LYS','881'), ('ILE','882'), ('VAL','883'), ('CYS','884'), ('GLN','885'),
            ('VAL','886'), ('GLY','887'), ('ARG','888'), ('LEU','889'), ('ASP','890'),
            ('ARG','891'), ('GLY','892'), ('LYS','893'), ('SER','894'), ('ALA','895'),
            ('ILE','896'), ('LEU','897'), ('TYR','898'), ('VAL','899'), ('LYS','900'),
            ('SER','901'), ('LEU','902'), ('LEU','903'), ('TRP','904'), ('THR','905'),
            ('GLN','906'), ('THR','907'), ('PHE','908'), ('MET','909'), ('ASN','910'),
            ('LYS','911'), ('GLU','912'), ('ASN','913'), ('GLN','914'), ('ASN','915'),
            ('HIS','916'), ('SER','917'), ('TYR','918'), ('SER','919'), ('LEU','920'),
            ('LYS','921'), ('SER','922'), ('SER','923'), ('ALA','924'), ('SER','925'),
            ('PHE','926'), ('ASN','927'), ('VAL','928'), ('ILE','929'), ('GLU','930'),
            ('PHE','931'), ('PRO','932'), ('TYR','933'), ('LYS','934'), ('ASN','935'),
            ('LEU','936'), ('PRO','937'), ('ILE','938'), ('GLU','939'), ('ASP','940'),
            ('ILE','941'), ('THR','942'), ('ASN','943'), ('SER','944'), ('THR','945'),
            ('LEU','946'), ('VAL','947'), ('THR','948'), ('THR','949'), ('ASN','950'),
            ('VAL','951'), ('THR','952'), ('TRP','953'), ('GLY','954'), ('ILE','955'),
            ('GLN','956'), ('PRO','957'), ('ALA','958'), ('PRO','959'), ('MET','960'),
            ('PRO','961'), ('VAL','962'), ('PRO','963'), ('VAL','964'), ('TRP','965'),
            ('VAL','966'), ('ILE','967')
            ] 

        self.xyzSequenceB = [
            ('GLY','1'), ('PRO','2'), ('ASN','3'),
            ('ILE','4'), ('CYS','5'), ('THR','6'), ('THR','7'), ('ARG','8'), ('GLY','9'),
            ('VAL','10'), ('SER','11'), ('SER','12'), ('CYS','13'), ('GLN','14'),
            ('GLN','15'), ('CYS','16'), ('LEU','17'), ('ALA','18'), ('VAL','19'),
            ('SER','20'), ('PRO','21'), ('MET','22'), ('CYS','23'), ('ALA','24'),
            ('TRP','25'), ('CYS','26'), ('SER','27'), ('ASP','28'), ('GLU','29'),
            ('ALA','30'), ('LEU','31'), ('PRO','32'), ('LEU','33'), ('GLY','34'),
            ('SER','35'), ('PRO','36'), ('ARG','37'), ('CYS','38'), ('ASP','39'),
            ('LEU','40'), ('LYS','41'), ('GLU','42'), ('ASN','43'), ('LEU','44'),
            ('LEU','45'), ('LYS','46'), ('ASP','47'), ('ASN','48'), ('CYS','49'),
            ('ALA','50'), ('PRO','51'), ('GLU','52'), ('SER','53'), ('ILE','54'),
            ('GLU','55'), ('PHE','56'), ('PRO','57'), ('VAL','58'), ('SER','59'),
            ('GLU','60'), ('ALA','61'), ('ARG','62'), ('VAL','63'), ('LEU','64'),
            ('GLU','65'), ('ASP','66'), ('ARG','67'), ('PRO','68'), ('LEU','69'),
            ('SER','70'), ('ASP','71'), ('LYS','72'), ('GLY','73'), ('SER','74'),
            ('GLY','75'), ('ASP','76'), ('SER','77'), ('SER','78'), ('GLN','79'),
            ('VAL','80'), ('THR','81'), ('GLN','82'), ('VAL','83'), ('SER','84'),
            ('PRO','85'), ('GLN','86'), ('ARG','87'), ('ILE','88'), ('ALA','89'),
            ('LEU','90'), ('ARG','91'), ('LEU','92'), ('ARG','93'), ('PRO','94'),
            ('ASP','95'), ('ASP','96'), ('SER','97'), ('LYS','98'), ('ASN','99'),
            ('PHE','100'), ('SER','101'), ('ILE','102'), ('GLN','103'), ('VAL','104'),
            ('ARG','105'), ('GLN','106'), ('VAL','107'), ('GLU','108'), ('ASP','109'),
            ('TYR','110'), ('PRO','111'), ('VAL','112'), ('ASP','113'), ('ILE','114'),
            ('TYR','115'), ('TYR','116'), ('LEU','117'), ('MET','118'), ('ASP','119'),
            ('LEU','120'), ('SER','121'), ('TYR','122'), ('SER','123'), ('MET','124'),
            ('LYS','125'), ('ASP','126'), ('ASP','127'), ('LEU','128'), ('TRP','129'),
            ('SER','130'), ('ILE','131'), ('GLN','132'), ('ASN','133'), ('LEU','134'),
            ('GLY','135'), ('THR','136'), ('LYS','137'), ('LEU','138'), ('ALA','139'),
            ('THR','140'), ('GLN','141'), ('MET','142'), ('ARG','143'), ('LYS','144'),
            ('LEU','145'), ('THR','146'), ('SER','147'), ('ASN','148'), ('LEU','149'),
            ('ARG','150'), ('ILE','151'), ('GLY','152'), ('PHE','153'), ('GLY','154'),
            ('ALA','155'), ('PHE','156'), ('VAL','157'), ('ASP','158'), ('LYS','159'),
            ('PRO','160'), ('VAL','161'), ('SER','162'), ('PRO','163'), ('TYR','164'),
            ('MET','165'), ('TYR','166'), ('ILE','167'), ('SER','168'), ('PRO','169'),
            ('PRO','170'), ('GLU','171'), ('ALA','172'), ('LEU','173'), ('GLU','174'),
            ('ASN','175'), ('PRO','176'), ('CYS','177'), ('TYR','178'), ('ASP','179'),
            ('MET','180'), ('LYS','181'), ('THR','182'), ('THR','183'), ('CYS','184'),
            ('LEU','185'), ('PRO','186'), ('MET','187'), ('PHE','188'), ('GLY','189'),
            ('TYR','190'), ('LYS','191'), ('HIS','192'), ('VAL','193'), ('LEU','194'),
            ('THR','195'), ('LEU','196'), ('THR','197'), ('ASP','198'), ('GLN','199'),
            ('VAL','200'), ('THR','201'), ('ARG','202'), ('PHE','203'), ('ASN','204'),
            ('GLU','205'), ('GLU','206'), ('VAL','207'), ('LYS','208'), ('LYS','209'),
            ('GLN','210'), ('SER','211'), ('VAL','212'), ('SER','213'), ('ARG','214'),
            ('ASN','215'), ('ARG','216'), ('ASP','217'), ('ALA','218'), ('PRO','219'),
            ('GLU','220'), ('GLY','221'), ('GLY','222'), ('PHE','223'), ('ASP','224'),
            ('ALA','225'), ('ILE','226'), ('MET','227'), ('GLN','228'), ('ALA','229'),
            ('THR','230'), ('VAL','231'), ('CYS','232'), ('ASP','233'), ('GLU','234'),
            ('LYS','235'), ('ILE','236'), ('GLY','237'), ('TRP','238'), ('ARG','239'),
            ('ASN','240'), ('ASP','241'), ('ALA','242'), ('SER','243'), ('HIS','244'),
            ('LEU','245'), ('LEU','246'), ('VAL','247'), ('PHE','248'), ('THR','249'),
            ('THR','250'), ('ASP','251'), ('ALA','252'), ('LYS','253'), ('THR','254'),
            ('HIS','255'), ('ILE','256'), ('ALA','257'), ('LEU','258'), ('ASP','259'),
            ('GLY','260'), ('ARG','261'), ('LEU','262'), ('ALA','263'), ('GLY','264'),
            ('ILE','265'), ('VAL','266'), ('GLN','267'), ('PRO','268'), ('ASN','269'),
            ('ASP','270'), ('GLY','271'), ('GLN','272'), ('CYS','273'), ('HIS','274'),
            ('VAL','275'), ('GLY','276'), ('SER','277'), ('ASP','278'), ('ASN','279'),
            ('HIS','280'), ('TYR','281'), ('SER','282'), ('ALA','283'), ('SER','284'),
            ('THR','285'), ('THR','286'), ('MET','287'), ('ASP','288'), ('TYR','289'),
            ('PRO','290'), ('SER','291'), ('LEU','292'), ('GLY','293'), ('LEU','294'),
            ('MET','295'), ('THR','296'), ('GLU','297'), ('LYS','298'), ('LEU','299'),
            ('SER','300'), ('GLN','301'), ('LYS','302'), ('ASN','303'), ('ILE','304'),
            ('ASN','305'), ('LEU','306'), ('ILE','307'), ('PHE','308'), ('ALA','309'),
            ('VAL','310'), ('THR','311'), ('GLU','312'), ('ASN','313'), ('VAL','314'),
            ('VAL','315'), ('ASN','316'), ('LEU','317'), ('TYR','318'), ('GLN','319'),
            ('ASN','320'), ('TYR','321'), ('SER','322'), ('GLU','323'), ('LEU','324'),
            ('ILE','325'), ('PRO','326'), ('GLY','327'), ('THR','328'), ('THR','329'),
            ('VAL','330'), ('GLY','331'), ('VAL','332'), ('LEU','333'), ('SER','334'),
            ('MET','335'), ('ASP','336'), ('SER','337'), ('SER','338'), ('ASN','339'),
            ('VAL','340'), ('LEU','341'), ('GLN','342'), ('LEU','343'), ('ILE','344'),
            ('VAL','345'), ('ASP','346'), ('ALA','347'), ('TYR','348'), ('GLY','349'),
            ('LYS','350'), ('ILE','351'), ('ARG','352'), ('SER','353'), ('LYS','354'),
            ('VAL','355'), ('GLU','356'), ('LEU','357'), ('GLU','358'), ('VAL','359'),
            ('ARG','360'), ('ASP','361'), ('LEU','362'), ('PRO','363'), ('GLU','364'),
            ('GLU','365'), ('LEU','366'), ('SER','367'), ('LEU','368'), ('SER','369'),
            ('PHE','370'), ('ASN','371'), ('ALA','372'), ('THR','373'), ('CYS','374'),
            ('LEU','375'), ('ASN','376'), ('ASN','377'), ('GLU','378'), ('VAL','379'),
            ('ILE','380'), ('PRO','381'), ('GLY','382'), ('LEU','383'), ('LYS','384'),
            ('SER','385'), ('CYS','386'), ('MET','387'), ('GLY','388'), ('LEU','389'),
            ('LYS','390'), ('ILE','391'), ('GLY','392'), ('ASP','393'), ('THR','394'),
            ('VAL','395'), ('SER','396'), ('PHE','397'), ('SER','398'), ('ILE','399'),
            ('GLU','400'), ('ALA','401'), ('LYS','402'), ('VAL','403'), ('ARG','404'),
            ('GLY','405'), ('CYS','406'), ('PRO','407'), ('GLN','408'), ('GLU','409'),
            ('LYS','410'), ('GLU','411'), ('LYS','412'), ('SER','413'), ('PHE','414'),
            ('THR','415'), ('ILE','416'), ('LYS','417'), ('PRO','418'), ('VAL','419'),
            ('GLY','420'), ('PHE','421'), ('LYS','422'), ('ASP','423'), ('SER','424'),
            ('LEU','425'), ('ILE','426'), ('VAL','427'), ('GLN','428'), ('VAL','429'),
            ('THR','430'), ('PHE','431'), ('ASP','432'), ('CYS','433'), ('ASP','434'),
            ('CYS','435'), ('ALA','436'), ('CYS','437'), ('GLN','438'), ('ALA','439'),
            ('GLN','440'), ('ALA','441'), ('GLU','442'), ('PRO','443'), ('ASN','444'),
            ('SER','445'), ('HIS','446'), ('ARG','447'), ('CYS','448'), ('ASN','449'),
            ('ASN','450'), ('GLY','451'), ('ASN','452'), ('GLY','453'), ('THR','454'),
            ('PHE','455'), ('GLU','456'), ('CYS','457'), ('GLY','458'), ('VAL','459'),
            ('CYS','460'), ('ARG','461'), ('CYS','462'), ('GLY','463'), ('PRO','464'),
            ('GLY','465'), ('TRP','466'), ('LEU','467'), ('GLY','468'), ('SER','469'),
            ('GLN','470'), ('CYS','471'), ('GLU','472'), ('CYS','473'), ('SER','474'),
            ('GLU','475'), ('GLU','476'), ('ASP','477'), ('TYR','478'), ('ARG','479'),
            ('PRO','480'), ('SER','481'), ('GLN','482'), ('GLN','483'), ('ASP','484'),
            ('GLU','485'), ('CYS','486'), ('SER','487'), ('PRO','488'), ('ARG','489'),
            ('GLU','490'), ('GLY','491'), ('GLN','492'), ('PRO','493'), ('VAL','494'),
            ('CYS','495'), ('SER','496'), ('GLN','497'), ('ARG','498'), ('GLY','499'),
            ('GLU','500'), ('CYS','501'), ('LEU','502'), ('CYS','503'), ('GLY','504'),
            ('GLN','505'), ('CYS','506'), ('VAL','507'), ('CYS','508'), ('HIS','509'),
            ('SER','510'), ('SER','511'), ('ASP','512'), ('PHE','513'), ('GLY','514'),
            ('LYS','515'), ('ILE','516'), ('THR','517'), ('GLY','518'), ('LYS','519'),
            ('TYR','520'), ('CYS','521'), ('GLU','522'), ('CYS','523'), ('ASP','524'),
            ('ASP','525'), ('PHE','526'), ('SER','527'), ('CYS','528'), ('VAL','529'),
            ('ARG','530'), ('TYR','531'), ('LYS','532'), ('GLY','533'), ('GLU','534'),
            ('MET','535'), ('CYS','536'), ('SER','537'), ('GLY','538'), ('HIS','539'),
            ('GLY','540'), ('GLN','541'), ('CYS','542'), ('SER','543'), ('CYS','544'),
            ('GLY','545'), ('ASP','546'), ('CYS','547'), ('LEU','548'), ('CYS','549'),
            ('ASP','550'), ('SER','551'), ('ASP','552'), ('TRP','553'), ('THR','554'),
            ('GLY','555'), ('TYR','556'), ('TYR','557'), ('CYS','558'), ('ASN','559'),
            ('CYS','560'), ('THR','561'), ('THR','562'), ('ARG','563'), ('THR','564'),
            ('ASP','565'), ('THR','566'), ('CYS','567'), ('MET','568'), ('SER','569'),
            ('SER','570'), ('ASN','571'), ('GLY','572'), ('LEU','573'), ('LEU','574'),
            ('CYS','575'), ('SER','576'), ('GLY','577'), ('ARG','578'), ('GLY','579'),
            ('LYS','580'), ('CYS','581'), ('GLU','582'), ('CYS','583'), ('GLY','584'),
            ('SER','585'), ('CYS','586'), ('VAL','587'), ('CYS','588'), ('ILE','589'),
            ('GLN','590'), ('PRO','591'), ('GLY','592'), ('SER','593'), ('TYR','594'),
            ('GLY','595'), ('ASP','596'), ('THR','597'), ('CYS','598'), ('GLU','599'),
            ('LYS','600'), ('CYS','601'), ('PRO','602'), ('THR','603'), ('CYS','604'),
            ('PRO','605'), ('ASP','606'), ('ALA','607'), ('CYS','608'), ('THR','609'),
            ('PHE','610'), ('LYS','611'), ('LYS','612'), ('GLU','613'), ('CYS','614'),
            ('VAL','615'), ('GLU','616'), ('CYS','617'), ('LYS','618'), ('LYS','619'),
            ('PHE','620'), ('ASP','621'), ('ARG','622'), ('GLY','623'), ('ALA','624'),
            ('LEU','625'), ('HIS','626'), ('ASP','627'), ('GLU','628'), ('ASN','629'),
            ('THR','630'), ('CYS','631'), ('ASN','632'), ('ARG','633'), ('TYR','634'),
            ('CYS','635'), ('ARG','636'), ('ASP','637'), ('GLU','638'), ('ILE','639'),
            ('GLU','640'), ('SER','641'), ('VAL','642'), ('LYS','643'), ('GLU','644'),
            ('LEU','645'), ('LYS','646'), ('ASP','647'), ('THR','648'), ('GLY','649'),
            ('LYS','650'), ('ASP','651'), ('ALA','652'), ('VAL','653'), ('ASN','654'),
            ('CYS','655'), ('THR','656'), ('TYR','657'), ('LYS','658'), ('ASN','659'),
            ('GLU','660'), ('ASP','661'), ('ASP','662'), ('CYS','663'), ('VAL','664'),
            ('VAL','665'), ('ARG','666'), ('PHE','667'), ('GLN','668'), ('TYR','669'),
            ('TYR','670'), ('GLU','671'), ('ASP','672'), ('SER','673'), ('SER','674'),
            ('GLY','675'), ('LYS','676'), ('SER','677'), ('ILE','678'), ('LEU','679'),
            ('TYR','680'), ('VAL','681'), ('VAL','682'), ('GLU','683'), ('GLU','684'),
            ('PRO','685'), ('GLU','686'), ('CYS','687'), ('PRO','688'), ('LYS','689'),
            ('GLY','690'), ('PRO','691'), ('ASP','692'), ('ILE','693'), ('LEU','694'),
            ('VAL','695')
            ]

    def getIdCode(self):
        return self.idCode

    def getSeqIdList(self):
        return self.seqIdList

    def getAuthSequenceList(self,chainId):
        if chainId == "A":
            return self.authSequenceA.split()
        elif chainId == "B":
            return self.authSequenceB.split()
        else:
            return []

    def getAuthSequenceListTest(self,chainId):
        seqT=self.getAuthSequenceList(chainId)
        self.__mutateSequence(seqT)
        return seqT

    def getAuthSequenceTestWithIndexList(self,chainId):
        sL0=self.getAuthSequenceListTest(chainId)
        sL=[]
        ir = 1
        for aa in sL0:
            sL.append( (aa,str(ir),'',ir) )
            ir += 1                
        return sL
    
    def __mutateSequence(self,seqIn):
        if len(seqIn) < 1: return        
        try:
            lenS=len(seqIn)
            maxDelete=5
            maxInsert=5
            maxInsertSeq=2            
            # random insert
            iLen=random.randint(1,maxInsert)
            for i in range(1,iLen):
                idx=random.randint(0,lenS)
                aa=random.choice(self.aaList)
                self.__insert(seqIn,idx,aa)

            iLen=random.randint(1,maxInsertSeq)
            for i in range(1,iLen):
                idx=random.randint(0,lenS)
                aa=random.choice(self.aaList)
                for ii in range(1,random.randint(1,10)):
                    self.__insert(seqIn,idx,aa)

            #deletions < 5
            dLen=random.randint(1,maxDelete)
            for i in range(1,dLen):
                idx=random.randint(0,lenS)                
                self.__remove(seqIn,idx)
            
        except:
            pass

        return seqIn

    def __deletionMutationSequence(self,seqIn):
        if len(seqIn) < 1: return
        try:
            lenS=len(seqIn)
            maxDelete=5
            #deletions < 5
            dLen=random.randint(1,maxDelete)
            for i in range(1,dLen):
                idx=random.randint(0,lenS)                
                self.__remove(seqIn,idx)
        except:
            pass
        return seqIn
    
    def __remove(self,seq3List, index):
        try:
            seq3List.remove(index)
        except:
            pass

    def __replace(self,seq3List, index, sInp):
        try:
            seq3List[index]=sInp
        except:
            pass
        
    def __insert(self,seq3List, index, sInp):
        try:
            rL=sInp.split()
            rL.reverse()
            for r in rL:
                seq3List.insert(index,r)
        except:
            pass
        
        
    def getAuthSequenceWithIndexList(self,chainId):
        sL=[]
        ir = 1
        if chainId == "A":
            for aa in self.authSequenceA.split():
                sL.append( (aa,str(ir),'',ir) )
                ir += 1                
        elif chainId == "B":
            for aa in self.authSequenceB.split():
                sL.append( (aa,str(ir),'',ir) )
                ir += 1                
        else:
            pass

        return sL
    
    def getXyzSequenceWithIndexList(self,chainId):
        if chainId == "A":
            seqT=self.xyzSequenceA
        elif chainId == "B":
            seqT=self.xyzSequenceB
        else:
            seqT=[]
        oSeq=[]
        ir=1
        for tup in seqT:
            oSeq.append((tup[0],tup[1],'',ir))
            ir+=1
        return oSeq
        

    def getXyzSequenceTestWithIndexList(self,chainId):
        if chainId == "A":
            seqT=self.xyzSequenceAMod
        elif chainId == "B":
            seqT=self.xyzSequenceB
        else:
            seqT=[]
        seqT=self.__deletionMutationSequence(seqT)
        oSeq=[]
        ir=1
        for tup in seqT:
            oSeq.append((tup[0],tup[1],'',ir))
            ir+=1
        return oSeq
        
        return oSeq
    
    def toList(self,strIn):
        sL=[]
        for ss in strIn:
            if ss in string.whitespace:
                continue
            sL.append(ss)
        return sL


    def getRefFeatureDict(self,chainId):
        seqFeature=SequenceFeature()

        if chainId == "A":
            seqFeature.setId(dbName='UNIPROT',dbCode="ITAV_HUMAN", dbAccession="P06756")
            seqFeature.setSource(organism="Homo sapiens",strain='', taxid='')
            seqFeature.setItem('MATCH_LENGTH',len(self.refSequenceA.strip()))
            return seqFeature.get()
        elif chainId == "B":
            seqFeature.setId(dbName='UNIPROT',dbCode="ITB3_HUMAN", dbAccession="P05106")
            seqFeature.setSource(organism="Homo sapiens",strain='', taxid='')
            seqFeature.setItem('MATCH_LENGTH',len(self.refSequenceB.strip()))            
            return seqFeature.get()            
        else:
            seqFeature.setId(dbName='UNIPROT',dbCode="HUMAN_"+chainId, dbAccession="TR000"+chainId)
            seqFeature.setSource(organism="Homo sapiens",strain='', taxid='')
            seqFeature.setItem('MATCH_LENGTH',len(self.refSequenceA.strip()))            
            return seqFeature.get()            



    def getAuthFeatureDict(self,chainId):
        seqFeature=SequenceFeature()
        seqFeature.setId(dbName='PDB',dbCode='3IJE', dbAccession='3IJE')
        seqFeature.setSource(organism="Homo sapiens",strain='', taxid='')        
        return seqFeature.get()

    def getXyzFeatureDict(self,chainId):
        seqFeature=SequenceFeature()
        seqFeature.setId(dbName='PDB',dbCode='3IJE', dbAccession='3IJE')
        seqFeature.setSource(organism="Homo sapiens",strain='', taxid='')        
        return seqFeature.get()        

    
    def getRefSequenceList(self,chainId):
        if chainId == "A":
            return self.toList(self.refSequenceA.strip())
        elif chainId == "B":
            return self.toList(self.refSequenceB.strip())
        else:
            return []

    def getRefSequenceWithIndexList(self,chainId):
        sL=[]
        ir=0
        if chainId == "A":
            for aa in self.__srd.cnv1To3List(self.toList(self.refSequenceA.strip()),'AA'):
                sL.append( (aa,str(ir),''))
                ir += 1
        elif chainId == "B":
            for aa in self.__srd.cnv1To3List(self.toList(self.refSequenceB.strip()),'AA'):
                sL.append( (aa,str(ir),'') )
                ir += 1
        else:
            pass

        return sL

    def getRefSequenceTestWithIndexList(self,chainId):
        if chainId == "A":
            seqT=self.__srd.cnv1To3List(self.toList(self.refSequenceA.strip()),'AA')
        elif chainId == "B":
            seqT=self.__srd.cnv1To3List(self.toList(self.refSequenceB.strip()),'AA')
        else:
            seqT=[]
        seqT=self.__mutateSequence(seqT)
        #
        sL=[]
        ir=0
        for aa in seqT:
            sL.append( (aa,str(ir),''))
            ir += 1

        return sL
    

    def getRefSequence3List(self,chainId):
        return self.__srd.cnv1To3List(self.getRefSequenceList(chainId),'AA')

    

def formatCpp(sL,seqName):
    sLen=len(sL)
    ibeg=0
    iend=0
    nW=10
    sys.stdout.write("// %s Length = %s \n" % (seqName,sLen))
    sys.stdout.write("static const char *_%sSequence[] = {\n" % seqName)                    
    while (ibeg < sLen):
        iend = min(ibeg+nW,sLen)
        for r in sL[ibeg:iend]:
            sys.stdout.write('"%s",' % r)
        sys.stdout.write('\n')
        ibeg = iend
    sys.stdout.write('""};\n')            
    
if __name__ == '__main__':
    sE=SequenceExamples()
    #sL=sE.getRefSequenceWithIndexList('A')
    sL= sE.getRefSequence3List('A')
    formatCpp(sL,"refA")    
    sL=sE.getAuthSequenceList('A')
    formatCpp(sL,"authA")

    #
    # Test sequence with random insertions and deletions
    #
    
    sTests={}
    for tt in ['authAT1','authAT2','authAT3','authAT4','authAT5','authAT6','authAT7','authAT8','authAT9','authAT10' ]:
        sL=sE.getAuthSequenceListTest('A')
        formatCpp(sL,tt)
        

