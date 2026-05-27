import os
import tempfile
import importlib.util


def load_module():
    repo_root = os.path.dirname(os.path.dirname(__file__))
    mod_path = os.path.join(repo_root, 'scripts', 'create_wgs_calling_intervals.py')
    spec = importlib.util.spec_from_file_location('create_wgs_calling_intervals', mod_path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def write_nbed(content):
    f = tempfile.NamedTemporaryFile('w', delete=False)
    f.write(content)
    f.close()
    return f.name


def test_get_nonN_region_splitting():
    mod = load_module()
    fa = {'chr1': 100, 'chr2': 80}
    nbed = write_nbed('\n'.join([
        'chr1\t1\t10',
        'chr1\t50\t60',
        'chr2\t20\t30',
    ]) + '\n')
    try:
        reg = mod.get_nonN_region(nbed, 5, fa)
        assert reg['chr1'] == [[11, 49], [61, 100]]
        assert reg['chr2'] == [[1, 19], [31, 80]]
    finally:
        os.unlink(nbed)


def test_get_nonN_region_merge_small_N():
    mod = load_module()
    fa = {'chr1': 100, 'chr2': 80}
    nbed = write_nbed('\n'.join([
        'chr1\t1\t10',
        'chr1\t50\t60',
        'chr2\t20\t30',
    ]) + '\n')
    try:
        # With a very large threshold, treat all N as 'small' (ignored)
        reg = mod.get_nonN_region(nbed, 1000, fa)
        assert reg['chr1'] == [[1, 100]]
        assert reg['chr2'] == [[1, 80]]
    finally:
        os.unlink(nbed)
