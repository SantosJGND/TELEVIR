"""
Host genome library for TELE-Vir.

Provides access to host genome FTP configurations for host depletion.
Data is centrally managed in sources.yaml.

Usage:
    from host_library import Host, HomoSapiens
    
    # Class-based access (backward compatible)
    human = HomoSapiens()
    print(human.remote_host, human.remote_path)
    
    # Or use centralized loader directly
    from load_sources import get_host_config
    config = get_host_config('homo_sapiens')
"""

from typing import Optional, Dict
from install_scripts.load_sources import get_loader, get_host_config


class Host:
    """Base class for host genome configurations."""
    remote_host: str
    remote_path: str
    remote_filename: str
    host_name: str
    common_name: str

    @classmethod
    def from_config(cls, config: Dict[str, str]):
        """Create Host instance from config dict."""
        instance = cls()
        instance.remote_host = config.get('host', '')
        instance.remote_path = config.get('path', '')
        instance.remote_filename = config.get('file', '')
        instance.host_name = config.get('host_name', '')
        instance.common_name = config.get('common_name', '')
        return instance


class HomoSapiens(Host):
    def __init__(self, config: Optional[Dict] = None):
        if config is None:
            config = get_host_config('homo_sapiens')
        if config:
            self.remote_host = config.get('host', 'ftp.ncbi.nlm.nih.gov')
            self.remote_path = config.get('path', '/genomes/refseq/vertebrate_mammalian/Homo_sapiens/latest_assembly_versions/GCF_000001405.40_GRCh38.p14/')
            self.remote_filename = config.get('file', 'GCF_000001405.40_GRCh38.p14_genomic.fna.gz')
            self.host_name = config.get('host_name', 'hg38')
            self.common_name = config.get('common_name', 'human')
        else:
            self.remote_host = "ftp.ncbi.nlm.nih.gov"
            self.remote_path = "/genomes/refseq/vertebrate_mammalian/Homo_sapiens/latest_assembly_versions/GCF_000001405.40_GRCh38.p14/"
            self.remote_filename = "GCF_000001405.40_GRCh38.p14_genomic.fna.gz"
            self.host_name = "hg38"
            self.common_name = "human"


class MarmotaMarmota(Host):
    def __init__(self, config: Optional[Dict] = None):
        if config is None:
            config = get_host_config('marmota_marmota')
        if config:
            self.remote_host = config.get('host', 'ftp.ncbi.nlm.nih.gov')
            self.remote_path = config.get('path', 'genomes/refseq/vertebrate_mammalian/Marmota_marmota/latest_assembly_versions/GCF_001458135.2_marMar/')
            self.remote_filename = config.get('file', 'GCF_001458135.2_marMar_genomic.fna.gz')
            self.host_name = config.get('host_name', 'marmota_marmota')
            self.common_name = config.get('common_name', 'marmot')
        else:
            self.remote_host = "ftp.ncbi.nlm.nih.gov"
            self.remote_path = "genomes/refseq/vertebrate_mammalian/Marmota_marmota/latest_assembly_versions/GCF_001458135.2_marMar/"
            self.remote_filename = "GCF_001458135.2_marMar_genomic.fna.gz"
            self.host_name = "marmota_marmota"
            self.common_name = "marmot"


class CulexPipiens(Host):
    def __init__(self, config: Optional[Dict] = None):
        if config is None:
            config = get_host_config('culex_pipiens')
        if config:
            self.remote_host = config.get('host', 'ftp.ncbi.nlm.nih.gov')
            self.remote_path = config.get('path', 'genomes/refseq/invertebrate/Culex_pipiens/latest_assembly_versions/GCF_016801865.2_TS_CPP_V2/')
            self.remote_filename = config.get('file', 'GCF_016801865.2_TS_CPP_V2_genomic.fna.gz')
            self.host_name = config.get('host_name', 'culex_pipiens')
            self.common_name = config.get('common_name', 'mosquito')
        else:
            self.remote_host = "ftp.ncbi.nlm.nih.gov"
            self.remote_path = "genomes/refseq/invertebrate/Culex_pipiens/latest_assembly_versions/GCF_016801865.2_TS_CPP_V2/"
            self.remote_filename = "GCF_016801865.2_TS_CPP_V2_genomic.fna.gz"
            self.host_name = "culex_pipiens"
            self.common_name = "mosquito"


class FelisCatus(Host):
    def __init__(self, config: Optional[Dict] = None):
        if config is None:
            config = get_host_config('felis_catus')
        if config:
            self.remote_host = config.get('host', 'ftp.ncbi.nlm.nih.gov')
            self.remote_path = config.get('path', 'genomes/refseq/vertebrate_mammalian/Felis_catus/latest_assembly_versions/GCF_018350175.1_F.catus_Fca126_mat1.0/')
            self.remote_filename = config.get('file', 'GCF_018350175.1_F.catus_Fca126_mat1.0_genomic.fna.gz')
            self.host_name = config.get('host_name', 'felis_catus')
            self.common_name = config.get('common_name', 'cat')
        else:
            self.remote_host = "ftp.ncbi.nlm.nih.gov"
            self.remote_path = "genomes/refseq/vertebrate_mammalian/Felis_catus/latest_assembly_versions/GCF_018350175.1_F.catus_Fca126_mat1.0/"
            self.remote_filename = "GCF_018350175.1_F.catus_Fca126_mat1.0_genomic.fna.gz"
            self.host_name = "felis_catus"
            self.common_name = "cat"


class CanisLupusFamiliaris(Host):
    def __init__(self, config: Optional[Dict] = None):
        if config is None:
            config = get_host_config('canis_lupus_familiaris')
        if config:
            self.remote_host = config.get('host', 'ftp.ncbi.nlm.nih.gov')
            self.remote_path = config.get('path', 'genomes/refseq/vertebrate_mammalian/Canis_lupus_familiaris/latest_assembly_versions/GCF_011100685.1_UU_Cfam_GSD_1.0/')
            self.remote_filename = config.get('file', 'GCF_011100685.1_UU_Cfam_GSD_1.0_genomic.fna.gz')
            self.host_name = config.get('host_name', 'canis_lupus_familiaris')
            self.common_name = config.get('common_name', 'dog')
        else:
            self.remote_host = "ftp.ncbi.nlm.nih.gov"
            self.remote_path = "genomes/refseq/vertebrate_mammalian/Canis_lupus_familiaris/latest_assembly_versions/GCF_011100685.1_UU_Cfam_GSD_1.0/"
            self.remote_filename = "GCF_011100685.1_UU_Cfam_GSD_1.0_genomic.fna.gz"
            self.host_name = "canis_lupus_familiaris"
            self.common_name = "dog"


class CyprinusCarpio(Host):
    def __init__(self, config: Optional[Dict] = None):
        if config is None:
            config = get_host_config('cyprinus_carpio')
        if config:
            self.remote_host = config.get('host', 'ftp.ncbi.nlm.nih.gov')
            self.remote_path = config.get('path', 'genomes/refseq/vertebrate_other/Cyprinus_carpio/latest_assembly_versions/GCF_018340385.1_ASM1834038v1/')
            self.remote_filename = config.get('file', 'GCF_018340385.1_ASM1834038v1_genomic.fna.gz')
            self.host_name = config.get('host_name', 'cyprinus_carpio')
            self.common_name = config.get('common_name', 'carp')
        else:
            self.remote_host = "ftp.ncbi.nlm.nih.gov"
            self.remote_path = "genomes/refseq/vertebrate_other/Cyprinus_carpio/latest_assembly_versions/GCF_018340385.1_ASM1834038v1/"
            self.remote_filename = "GCF_018340385.1_ASM1834038v1_genomic.fna.gz"
            self.host_name = "cyprinus_carpio"
            self.common_name = "carp"


class AnasPlatyrhynchos(Host):
    def __init__(self, config: Optional[Dict] = None):
        if config is None:
            config = get_host_config('anas_platyrhynchos')
        if config:
            self.remote_host = config.get('host', 'ftp.ncbi.nlm.nih.gov')
            self.remote_path = config.get('path', 'genomes/refseq/vertebrate_other/Anas_platyrhynchos/latest_assembly_versions/GCF_015476345.1_ZJU1.0/')
            self.remote_filename = config.get('file', 'GCF_015476345.1_ZJU1.0_genomic.fna.gz')
            self.host_name = config.get('host_name', 'anas_platyrhynchos')
            self.common_name = config.get('common_name', 'duck')
        else:
            self.remote_host = "ftp.ncbi.nlm.nih.gov"
            self.remote_path = "genomes/refseq/vertebrate_other/Anas_platyrhynchos/latest_assembly_versions/GCF_015476345.1_ZJU1.0/"
            self.remote_filename = "GCF_015476345.1_ZJU1.0_genomic.fna.gz"
            self.host_name = "anas_platyrhynchos"
            self.common_name = "duck"


class PipistrellusKuhlii(Host):
    def __init__(self, config: Optional[Dict] = None):
        if config is None:
            config = get_host_config('pipistrellus_kuhlii')
        if config:
            self.remote_host = config.get('host', 'ftp.ncbi.nlm.nih.gov')
            self.remote_path = config.get('path', 'genomes/refseq/vertebrate_mammalian/Pipistrellus_kuhlii/latest_assembly_versions/GCF_014108245.1_mPipKuh1.p/')
            self.remote_filename = config.get('file', 'GCF_014108245.1_mPipKuh1.p_genomic.fna.gz')
            self.host_name = config.get('host_name', 'pipistrellus_kuhlii')
            self.common_name = config.get('common_name', 'bat')
        else:
            self.remote_host = "ftp.ncbi.nlm.nih.gov"
            self.remote_path = "genomes/refseq/vertebrate_mammalian/Pipistrellus_kuhlii/latest_assembly_versions/GCF_014108245.1_mPipKuh1.p/"
            self.remote_filename = "GCF_014108245.1_mPipKuh1.p_genomic.fna.gz"
            self.host_name = "pipistrellus_kuhlii"
            self.common_name = "bat"


class PhlebotomusPapatasi(Host):
    def __init__(self, config: Optional[Dict] = None):
        if config is None:
            config = get_host_config('phlebotomus_papatasi')
        if config:
            self.remote_host = config.get('host', 'ftp.ncbi.nlm.nih.gov')
            self.remote_path = config.get('path', 'genomes/refseq/invertebrate/Phlebotomus_papatasi/latest_assembly_versions/GCF_024763615.1_Ppap_2.1/')
            self.remote_filename = config.get('file', 'GCF_024763615.1_Ppap_2.1_genomic.fna.gz')
            self.host_name = config.get('host_name', 'phlebotomus_papatasi')
            self.common_name = config.get('common_name', 'sandfly')
        else:
            self.remote_host = "ftp.ncbi.nlm.nih.gov"
            self.remote_path = "genomes/refseq/invertebrate/Phlebotomus_papatasi/latest_assembly_versions/GCF_024763615.1_Ppap_2.1/"
            self.remote_filename = "GCF_024763615.1_Ppap_2.1_genomic.fna.gz"
            self.host_name = "phlebotomus_papatasi"
            self.common_name = "sandfly"


class SusScrofa(Host):
    def __init__(self, config: Optional[Dict] = None):
        if config is None:
            config = get_host_config('sus_scrofa')
        if config:
            self.remote_host = config.get('host', 'ftp.ncbi.nlm.nih.gov')
            self.remote_path = config.get('path', '/genomes/refseq/vertebrate_mammalian/Sus_scrofa/latest_assembly_versions/GCF_000003025.6_Sscrofa11.1/')
            self.remote_filename = config.get('file', 'GCF_000003025.6_Sscrofa11.1_genomic.fna.gz')
            self.host_name = config.get('host_name', 'sus_scrofa')
            self.common_name = config.get('common_name', 'pig')
        else:
            self.remote_host = "ftp.ncbi.nlm.nih.gov"
            self.remote_path = "/genomes/refseq/vertebrate_mammalian/Sus_scrofa/latest_assembly_versions/GCF_000003025.6_Sscrofa11.1/"
            self.remote_filename = "GCF_000003025.6_Sscrofa11.1_genomic.fna.gz"
            self.host_name = "sus_scrofa"
            self.common_name = "pig"


class AedesAlbopictus(Host):
    def __init__(self, config: Optional[Dict] = None):
        if config is None:
            config = get_host_config('aedes_albopictus')
        if config:
            self.remote_host = config.get('host', 'ftp.ncbi.nlm.nih.gov')
            self.remote_path = config.get('path', 'genomes/refseq/invertebrate/Aedes_albopictus/latest_assembly_versions/GCF_035046485.1_AalbF5/')
            self.remote_filename = config.get('file', 'GCF_035046485.1_AalbF5_genomic.fna.gz')
            self.host_name = config.get('host_name', 'aedes_albopictus')
            self.common_name = config.get('common_name', 'mosquito')
        else:
            self.remote_host = "ftp.ncbi.nlm.nih.gov"
            self.remote_path = "genomes/refseq/invertebrate/Aedes_albopictus/latest_assembly_versions/GCF_035046485.1_AalbF5/"
            self.remote_filename = "GCF_035046485.1_AalbF5_genomic.fna.gz"
            self.host_name = "aedes_albopictus"
            self.common_name = "mosquito"


class GallusGallus(Host):
    def __init__(self, config: Optional[Dict] = None):
        if config is None:
            config = get_host_config('gallus_gallus')
        if config:
            self.remote_host = config.get('host', 'ftp.ncbi.nlm.nih.gov')
            self.remote_path = config.get('path', 'genomes/refseq/vertebrate_other/Gallus_gallus/latest_assembly_versions/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b/')
            self.remote_filename = config.get('file', 'GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.fna.gz')
            self.host_name = config.get('host_name', 'gallus_gallus')
            self.common_name = config.get('common_name', 'chicken')
        else:
            self.remote_host = "ftp.ncbi.nlm.nih.gov"
            self.remote_path = "genomes/refseq/vertebrate_other/Gallus_gallus/latest_assembly_versions/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b/"
            self.remote_filename = "GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.fna.gz"
            self.host_name = "gallus_gallus"
            self.common_name = "chicken"


class OncorhynchusMykiss(Host):
    def __init__(self, config: Optional[Dict] = None):
        if config is None:
            config = get_host_config('oncorhynchus_mykiss')
        if config:
            self.remote_host = config.get('host', 'ftp.ncbi.nlm.nih.gov')
            self.remote_path = config.get('path', 'genomes/refseq/vertebrate_other/Oncorhynchus_mykiss/latest_assembly_versions/GCF_013265735.2_USDA_OmykA_1.1/')
            self.remote_filename = config.get('file', 'GCF_013265735.2_USDA_OmykA_1.1_genomic.fna.gz')
            self.host_name = config.get('host_name', 'oncorhynchus_mykiss')
            self.common_name = config.get('common_name', 'rainbow_trout')
        else:
            self.remote_host = "ftp.ncbi.nlm.nih.gov"
            self.remote_path = "genomes/refseq/vertebrate_other/Oncorhynchus_mykiss/latest_assembly_versions/GCF_013265735.2_USDA_OmykA_1.1/"
            self.remote_filename = "GCF_013265735.2_USDA_OmykA_1.1_genomic.fna.gz"
            self.host_name = "oncorhynchus_mykiss"
            self.common_name = "rainbow_trout"


class SalmoSalar(Host):
    def __init__(self, config: Optional[Dict] = None):
        if config is None:
            config = get_host_config('salmo_salar')
        if config:
            self.remote_host = config.get('host', 'ftp.ncbi.nlm.nih.gov')
            self.remote_path = config.get('path', 'genomes/refseq/vertebrate_other/Salmo_salar/latest_assembly_versions/GCF_905237065.1_Ssal_v3.1/')
            self.remote_filename = config.get('file', 'GCF_905237065.1_Ssal_v3.1_genomic.fna.gz')
            self.host_name = config.get('host_name', 'salmo_salar')
            self.common_name = config.get('common_name', 'atlantic_salmon')
        else:
            self.remote_host = "ftp.ncbi.nlm.nih.gov"
            self.remote_path = "genomes/refseq/vertebrate_other/Salmo_salar/latest_assembly_versions/GCF_905237065.1_Ssal_v3.1/"
            self.remote_filename = "GCF_905237065.1_Ssal_v3.1_genomic.fna.gz"
            self.host_name = "salmo_salar"
            self.common_name = "atlantic_salmon"


class BosTaurus(Host):
    def __init__(self, config: Optional[Dict] = None):
        if config is None:
            config = get_host_config('bos_taurus')
        if config:
            self.remote_host = config.get('host', 'ftp.ncbi.nlm.nih.gov')
            self.remote_path = config.get('path', 'genomes/refseq/vertebrate_mammalian/Bos_taurus/latest_assembly_versions/GCF_002263795.3_ARS-UCD2.0/')
            self.remote_filename = config.get('file', 'GCF_002263795.3_ARS-UCD2.0_genomic.fna.gz')
            self.host_name = config.get('host_name', 'bos_taurus')
            self.common_name = config.get('common_name', 'cow')
        else:
            self.remote_host = "ftp.ncbi.nlm.nih.gov"
            self.remote_path = "genomes/refseq/vertebrate_mammalian/Bos_taurus/latest_assembly_versions/GCF_002263795.3_ARS-UCD2.0/"
            self.remote_filename = "GCF_002263795.3_ARS-UCD2.0_genomic.fna.gz"
            self.host_name = "bos_taurus"
            self.common_name = "cow"


class NeogaleVison(Host):
    def __init__(self, config: Optional[Dict] = None):
        if config is None:
            config = get_host_config('neogale_vison')
        if config:
            self.remote_host = config.get('host', 'ftp.ncbi.nlm.nih.gov')
            self.remote_path = config.get('path', 'genomes/refseq/vertebrate_mammalian/Neogale_vison/latest_assembly_versions/GCF_020171115.1_ASM_NN_V1/')
            self.remote_filename = config.get('file', 'GCF_020171115.1_ASM_NN_V1_genomic.fna.gz')
            self.host_name = config.get('host_name', 'neogale_vison')
            self.common_name = config.get('common_name', 'mink')
        else:
            self.remote_host = "ftp.ncbi.nlm.nih.gov"
            self.remote_path = "genomes/refseq/vertebrate_mammalian/Neogale_vison/latest_assembly_versions/GCF_020171115.1_ASM_NN_V1/"
            self.remote_filename = "GCF_020171115.1_ASM_NN_V1_genomic.fna.gz"
            self.host_name = "neogale_vison"
            self.common_name = "mink"


HOST_CLASSES = {
    'hg38': HomoSapiens,
    'homo_sapiens': HomoSapiens,
    'marmota_marmota': MarmotaMarmota,
    'culex_pipiens': CulexPipiens,
    'felis_catus': FelisCatus,
    'canis_lupus_familiaris': CanisLupusFamiliaris,
    'cyprinus_carpio': CyprinusCarpio,
    'anas_platyrhynchos': AnasPlatyrhynchos,
    'pipistrellus_kuhlii': PipistrellusKuhlii,
    'phlebotomus_papatasi': PhlebotomusPapatasi,
    'sus_scrofa': SusScrofa,
    'aedes_albopictus': AedesAlbopictus,
    'gallus_gallus': GallusGallus,
    'oncorhynchus_mykiss': OncorhynchusMykiss,
    'salmo_salar': SalmoSalar,
    'bos_taurus': BosTaurus,
    'neogale_vison': NeogaleVison,
}


def get_host(host_name: str) -> Optional[Host]:
    """Get host instance by name (looks up by host_name or common_name)."""
    for cls in Host.__subclasses__():
        instance = cls()
        if instance.host_name == host_name or instance.common_name == host_name:
            return instance
    return None
