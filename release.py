"""Release and publish ilus to PyPI.

Author: Shujia Huang
Date: 2021-04-30
"""
import importlib
from subprocess import call

spec = importlib.util.spec_from_file_location("_", "./setup.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)

# Publish ilus to PyPI by using twine: https://segmentfault.com/a/1190000008663126
# pip install twine
# python setup.py register
# python setup.py sdist build && twine upload dist/ilus-0.0.5.tar.gz

#call(["pandoc", "--from=markdown", "--to=rst", "-o", "README.rst", "README.md"])
call(["python", "setup.py", "sdist"])
tarball = "dist/{}-{}.tar.gz".format(module.meta.__DISTNAME__, module.meta.__VERSION__)
call(["twine", "upload", "-r", "pypi", tarball])

