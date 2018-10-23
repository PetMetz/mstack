import os
from setuptools import setup, find_packages

def fread(filename):
  return open(os.path.join(os.path.dirname(__file__), filename)).read()

def get_version():
  return fread('VERSION')

setup_args = dict(
	name = 'mstack',
	version = get_version(),
	license = 'BSD-3-Clause',
	author = 'Peter C. Metz',
	author_email = 'pcm1@alfred.edu',
	description = 'Stacking disorder tools for Python.',
	long_description = fread('README.md'),
	install_requires = [
		'scipy>=0.18.1',
		'numpy>=1.12.0',
		'matplotlib >=2.0.0',
		'lmfit>=0.9.5',
		'emcee>=2.2.1',
		'corner>=2.0.1',
        'schwimmbad>=0.3.0',
        'tabulate>=0.8.2',
        'dill>=0.2.7.1'


						],
	classifiers = [
		'Programming Language :: Python :: 2.7',
		'Operating System :: POSIX',
		'License :: OSI Approved :: BSD-3-Clause'
				   ],
	platforms = 'linux',
	packages = ['mstack'],
	include_package_data = True
	# find_packages(),
	# package_dir={'': 'mstack'},
	# include_package_data = True
		        )


if __name__ == '__main__':
	setup(**setup_args)

# EOF #
