"""
Unit tests for host_library.py - Host genome classes.
"""

import pytest
from host_library import (
    Host, HomoSapiens, MarmotaMarmota, CulexPipiens, FelisCatus,
    CanisLupusFamiliaris, CyprinusCarpio, AnasPlatyrhynchos,
    PipistrellusKuhlii, PhlebotomusPapatasi, SusScrofa, AedesAlbopictus,
    GallusGallus, OncorhynchusMykiss, SalmoSalar, BosTaurus, NeogaleVison,
    get_host, HOST_CLASSES
)


class TestHostBaseClass:
    """Tests for Host base class."""

    def test_host_class_exists(self):
        """Test Host base class exists."""
        assert Host is not None

    def test_host_from_config(self, sample_host_config):
        """Test Host.from_config creates valid instance."""
        host = Host.from_config(sample_host_config)
        assert host.host_name == sample_host_config['host_name']
        assert host.common_name == sample_host_config['common_name']
        assert host.remote_host == sample_host_config['host']
        assert host.remote_path == sample_host_config['path']
        assert host.remote_filename == sample_host_config['file']


class TestHomoSapiens:
    """Tests for HomoSapiens host class."""

    def test_homo_sapiens_attributes(self):
        """Test HomoSapiens has correct attributes."""
        human = HomoSapiens()
        assert human.host_name == 'hg38'
        assert human.common_name == 'human'
        assert 'ftp.ncbi' in human.remote_host
        assert 'Homo_sapiens' in human.remote_path
        assert human.remote_filename.endswith('.fna.gz')

    def test_homo_sapiens_ftp_path(self):
        """Test HomoSapiens FTP path is valid."""
        human = HomoSapiens()
        expected_path = "/genomes/refseq/vertebrate_mammalian/Homo_sapiens/latest_assembly_versions/GCF_000001405.40_GRCh38.p14/"
        assert human.remote_path == expected_path


class TestAllHostClasses:
    """Tests for all host genome classes."""

    @pytest.mark.parametrize("host_class,name,expected_name,expected_common", [
        (HomoSapiens, 'HomoSapiens', 'hg38', 'human'),
        (MarmotaMarmota, 'MarmotaMarmota', 'marmota_marmota', 'marmot'),
        (CulexPipiens, 'CulexPipiens', 'culex_pipiens', 'mosquito'),
        (FelisCatus, 'FelisCatus', 'felis_catus', 'cat'),
        (CanisLupusFamiliaris, 'CanisLupusFamiliaris', 'canis_lupus_familiaris', 'dog'),
        (CyprinusCarpio, 'CyprinusCarpio', 'cyprinus_carpio', 'carp'),
        (AnasPlatyrhynchos, 'AnasPlatyrhynchos', 'anas_platyrhynchos', 'duck'),
        (PipistrellusKuhlii, 'PipistrellusKuhlii', 'pipistrellus_kuhlii', 'bat'),
        (PhlebotomusPapatasi, 'PhlebotomusPapatasi', 'phlebotomus_papatasi', 'sandfly'),
        (SusScrofa, 'SusScrofa', 'sus_scrofa', 'pig'),
        (AedesAlbopictus, 'AedesAlbopictus', 'aedes_albopictus', 'mosquito'),
        (GallusGallus, 'GallusGallus', 'gallus_gallus', 'chicken'),
        (OncorhynchusMykiss, 'OncorhynchusMykiss', 'oncorhynchus_mykiss', 'rainbow_trout'),
        (SalmoSalar, 'SalmoSalar', 'salmo_salar', 'atlantic_salmon'),
        (BosTaurus, 'BosTaurus', 'bos_taurus', 'cow'),
        (NeogaleVison, 'NeogaleVison', 'neogale_vison', 'mink'),
    ])
    def test_host_class_attributes(self, host_class, name, expected_name, expected_common):
        """Test each host class has correct attributes."""
        host = host_class()
        assert host.host_name == expected_name, f"{name}: host_name mismatch"
        assert host.common_name == expected_common, f"{name}: common_name mismatch"
        assert host.remote_host == 'ftp.ncbi.nlm.nih.gov'
        assert host.remote_filename.endswith('.fna.gz')

    def test_all_host_classes_instantiate(self):
        """Test all host classes can be instantiated without error."""
        subclasses = Host.__subclasses__()
        for cls in subclasses:
            instance = cls()
            assert instance.host_name is not None
            assert instance.common_name is not None
            assert instance.remote_host is not None
            assert instance.remote_path is not None
            assert instance.remote_filename is not None

    def test_host_class_count(self):
        """Test that all 16 expected host classes exist."""
        subclasses = Host.__subclasses__()
        assert len(subclasses) == 16, f"Expected 16 host classes, found {len(subclasses)}"


class TestGetHost:
    """Tests for get_host function."""

    def test_get_host_by_host_name(self):
        """Test get_host by host_name (e.g., 'hg38')."""
        host = get_host('hg38')
        assert host is not None
        assert host.host_name == 'hg38'
        assert host.common_name == 'human'

    def test_get_host_by_common_name(self):
        """Test get_host by common_name (e.g., 'human')."""
        host = get_host('human')
        assert host is not None
        assert host.common_name == 'human'
        assert host.host_name == 'hg38'

    def test_get_host_pig(self):
        """Test get_host for pig."""
        host = get_host('pig')
        assert host is not None
        assert host.common_name == 'pig'
        assert host.host_name == 'sus_scrofa'

    def test_get_host_mosquito(self):
        """Test get_host for mosquito (multiple hosts)."""
        host = get_host('mosquito')
        assert host is not None
        assert host.common_name == 'mosquito'

    def test_get_host_invalid(self):
        """Test get_host returns None for invalid host."""
        host = get_host('invalid_host_name_xyz')
        assert host is None

    def test_get_host_none(self):
        """Test get_host returns None for None input."""
        host = get_host(None)
        assert host is None


class TestHostClassMapping:
    """Tests for HOST_CLASSES mapping."""

    def test_host_classes_dict_exists(self):
        """Test HOST_CLASSES dict exists."""
        assert HOST_CLASSES is not None
        assert isinstance(HOST_CLASSES, dict)

    def test_host_classes_has_required_keys(self):
        """Test HOST_CLASSES has expected keys."""
        required_keys = ['hg38', 'homo_sapiens', 'sus_scrofa']
        for key in required_keys:
            assert key in HOST_CLASSES, f"Missing key: {key}"

    def test_host_classes_values_are_classes(self):
        """Test HOST_CLASSES values are Host subclasses."""
        for key, cls in HOST_CLASSES.items():
            assert issubclass(cls, Host), f"Value for {key} is not a Host subclass"


class TestHostAttributeValidation:
    """Tests for validating host attribute formats."""

    def test_all_hosts_use_ncbi_ftp(self):
        """Test all hosts use NCBI FTP server."""
        for cls in Host.__subclasses__():
            host = cls()
            assert 'ncbi' in host.remote_host.lower(), \
                f"{cls.__name__} doesn't use NCBI FTP"

    def test_all_host_files_are_gzipped(self):
        """Test all host genome files are gzip compressed."""
        for cls in Host.__subclasses__():
            host = cls()
            assert host.remote_filename.endswith('.fna.gz'), \
                f"{cls.__name__} file is not .fna.gz: {host.remote_filename}"

    def test_all_host_paths_contain_genomes(self):
        """Test all host paths contain 'genomes'."""
        for cls in Host.__subclasses__():
            host = cls()
            assert 'genomes' in host.remote_path.lower(), \
                f"{cls.__name__} path doesn't contain 'genomes'"

    def test_all_host_paths_start_correctly(self):
        """Test all host paths start with expected structure."""
        for cls in Host.__subclasses__():
            host = cls()
            path = host.remote_path.lstrip('/')
            assert path.startswith('genomes/'), \
                f"{cls.__name__} path doesn't start with 'genomes/'"


class TestBackwardCompatibility:
    """Tests to ensure backward compatibility with old Host class usage."""

    def test_instantiation_without_config(self):
        """Test host classes can be instantiated without config."""
        human = HomoSapiens()
        assert human.host_name is not None

    def test_subclass_check(self):
        """Test Host subclasses are properly registered."""
        subclasses = Host.__subclasses__()
        assert HomoSapiens in subclasses
        assert SusScrofa in subclasses

    def test_host_lookup_consistency(self):
        """Test get_host returns consistent results with class instantiation."""
        for cls in Host.__subclasses__():
            instance = cls()
            found = get_host(instance.host_name)
            assert found is not None
            assert found.host_name == instance.host_name
