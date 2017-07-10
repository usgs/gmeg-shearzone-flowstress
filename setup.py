try:
	from setuptools import setup
except ImportError:
	from distutils.core import setup

config = {
	'description':'Flow stress calculator',
	'author': 'Ben Melosh',
	'url':'https://github.com/BenjiDa/Flowstress',
	'author email':'bmelosh@usgs.gov',
	'version': '0.1',
	'install_requires':['NAME'],
	'name':'Flow stress calculator'
}

setup(**config)

