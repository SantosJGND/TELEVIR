HOSTS_TO_INSTALL = [
    "hg38",
    "sus_scrofa",
    # "aedes_albopictus",
    # "gallus_gallus",
    # "oncorhynchus_mykiss",
    # "salmo_salar",
    # "bos_taurus",
    # "neogale_vison",
]


class Host:
    remote_host: str
    remote_path: str
    remote_filename: str
    host_name: str
    common_name: str


class Hg38(Host):
    def __init__(self):
        self.remote_host = "hgdownload.soe.ucsc.edu"
        self.remote_path = "/goldenPath/hg38/bigZips/"
        self.remote_filename = "hg38.fa.gz"
        self.host_name = "hg38"
        self.common_name = "human"


class SusScrofa(Host):
    def __init__(self):
        self.remote_host = "ftp.ncbi.nlm.nih.gov"
        self.remote_path = "genomes/refseq/vertebrate_mammalian/Sus_scrofa/latest_assembly_versions/GCF_000003025.6_Sscrofa11.1/"
        self.remote_filename = "GCF_000003025.6_Sscrofa11.1_genomic.fna.gz"
        self.host_name = "sus_scrofa"
        self.common_name = "pig"


class AedesAlbopictus(Host):
    def __init__(self):
        self.remote_host = "ftp.ncbi.nlm.nih.gov"
        self.remote_path = "genomes/refseq/invertebrate/Aedes_albopictus/latest_assembly_versions/GCF_006496715.2_Aalbo_primary.1/"
        self.remote_filename = "GCF_006496715.2_Aalbo_primary.1_genomic.fna.gz"
        self.host_name = "aedes_albopictus"
        self.common_name = "mosquito"


class GallusGallus(Host):
    def __init__(self):
        self.remote_host = "ftp.ncbi.nlm.nih.gov"
        self.remote_path = "genomes/refseq/vertebrate_other/Gallus_gallus/latest_assembly_versions/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b/"
        self.remote_filename = (
            "GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.fna.gz"
        )
        self.host_name = "gallus_gallus"
        self.common_name = "chicken"


class OncorhynchusMykiss(Host):
    def __init__(self):
        self.remote_host = "ftp.ncbi.nlm.nih.gov"
        self.remote_path = "genomes/refseq/vertebrate_other/Oncorhynchus_mykiss/latest_assembly_versions/GCF_013265735.2_USDA_OmykA_1.1/"
        self.remote_filename = "GCF_013265735.2_USDA_OmykA_1.1_genomic.fna.gz"
        self.host_name = "oncorhynchus_mykiss"
        self.common_name = "rainbow_trout"


class SalmoSalar(Host):
    def __init__(self):
        self.remote_host = "ftp.ncbi.nlm.nih.gov"
        self.remote_path = "genomes/refseq/vertebrate_other/Acipenser_ruthenus/latest_assembly_versions/GCF_902713425.1_fAciRut3.2_maternal_haplotype/"
        self.remote_filename = (
            "GCF_902713425.1_fAciRut3.2_maternal_haplotype_genomic.fna.gz"
        )
        self.host_name = "salmo_salar"
        self.common_name = "atlantic_salmon"


class BosTaurus(Host):
    def __init__(self):
        self.remote_host = "ftp.ncbi.nlm.nih.gov"
        self.remote_path = "genomes/refseq/vertebrate_mammalian/Bos_taurus/latest_assembly_versions/GCF_002263795.2_ARS-UCD1.2/"
        self.remote_filename = "GCF_002263795.2_ARS-UCD1.2_genomic.fna.gz"
        self.host_name = "bos_taurus"
        self.common_name = "cow"


class NeogaleVison(Host):
    def __init__(self):
        self.remote_host = "ftp.ncbi.nlm.nih.gov"
        self.remote_path = "genomes/refseq/vertebrate_other/Neogale_vison/latest_assembly_versions/GCF_020171115.1_ASM_NN_V1/"
        self.remote_filename = "GCF_020171115.1_ASM_NN_V1_cds_from_genomic.fna.gz"
        self.host_name = "neogale_vison"
        self.common_name = "mink"
