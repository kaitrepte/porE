import setuptools
from numpy.distutils.core import Extension, setup

with open("README.md", "r") as fh:
     name="porE",
     long_description = fh.read()

setup(
   name="porE",
   version="1.0.4",
   author="Kai Trepte",
   author_email="kai.trepte1987@gmail.com",
   description="Porosity Evaluation tool",
   url="https://github.com/kaitrepte/porE", 
   license='APACHE2.0',
   long_description=long_description,
   long_description_content_type="text/markdown",
   include_package_data=True,
   packages = ['porE/lib','porE/gui','porE/io','porE/hea'],
   install_requires=['ase'],
   zip_safe=False,
   ext_modules=[Extension(name='pore', sources=['porE/lib/porE.f90'], f2py_options=['--quiet'])]
)
