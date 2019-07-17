import subprocess
import tempfile


def _exec_notebook(path):
    with tempfile.NamedTemporaryFile(suffix=".ipynb") as fout:
        args = ["jupyter", "nbconvert", "--to", "notebook", "--execute",
                "--ExecutePreprocessor.timeout=1000",
                "--ExecutePreprocessor.kernel_name=python3",
                "--output", fout.name, path]
        subprocess.check_call(args)


def test():
    print('Testing Jupyter notebooks...')
    _exec_notebook('QEr_2D_joint.ipynb')
    _exec_notebook('ERNR_bands.ipynb')
    _exec_notebook('ms_correction.ipynb')
    _exec_notebook('yield_width_compare.ipynb')
    _exec_notebook('fitting_errors.ipynb')
    _exec_notebook('EpEq_2D_joint.ipynb')
    _exec_notebook('EpEq_calc_contours.ipynb')
