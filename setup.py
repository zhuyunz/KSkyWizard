import setuptools

# Get some values from the setup.cfg
try:
    from ConfigParser import ConfigParser
except ImportError:
    from configparser import ConfigParser

conf = ConfigParser()
conf.read(['setup.cfg'])
metadata = dict(conf.items("metadata"))

NAME = metadata['name']
VERSION = metadata['version']
RELEASE = 'dev' not in VERSION
AUTHOR = metadata["author"]
AUTHOR_EMAIL = metadata["author_email"]
LICENSE = metadata["license"]
DESCRIPTION = metadata["description"]

entry_points = {
    'gui_scripts': [
        "kskywizard = kskywizard.kskywizard:main"
    ]}

setuptools.setup(name=NAME,
      provides=NAME,
      version=VERSION,
      license=LICENSE,
      description=DESCRIPTION,
      long_description=open('README.md').read(),
      author=AUTHOR,
      author_email=AUTHOR_EMAIL,
      packages=setuptools.find_packages(),
      package_data={'': ['data/extin/*', 'data/stds/*']},
      entry_points=entry_points,
      install_requires=[
        'astropy',
        'argparse',
        'matplotlib',
        'numpy',
        'PyQt5',
        'scipy']
)



