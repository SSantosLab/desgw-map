from setuptools import setup
import sandman

long_description = read('README')

setup(
    name='desgw-map',
    version="1.2",
    description='Compute DECam observations of LIGO counterparts',
    url="",
    author='James Annis'
    license="GPL",
    packages=['desgw-map']
    install_requires=[
          'healpy',
          'fslalib',
    ],
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Topic :: Scientific/Engineering :: Astronomy",
    ],

)
