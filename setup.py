from setuptools import setup
import re

def readme():
    with open('README.rst') as f:
        return f.read()

def version():
    v="NA"
    with open("barnaba/_version.py") as f:
        for line in f:
            if re.match("^__version__ *= *",line):
                v=re.sub('^__version__ *= *"',"",line.strip())
                v=re.sub('".*',"",v)
    return v

setup(name='barnaba',
      description='analyze nucleic acid 3D structures and MD trajectories',
      version=version(),
      long_description=readme(),
      classifiers=[
          'Development Status :: 5 - Production/Stable',
          'Intended Audience :: Science/Research',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
          'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
          'Programming Language :: Python :: 3.6',
          'Programming Language :: Python :: 3.7',
          'Programming Language :: Python :: 3.8',
          'Programming Language :: Python :: 3.9',
          'Programming Language :: Python :: 3.10',
      ],      
      url='https://github.com/srnas/barnaba',
      author='Sandro Bottaro',
      author_email='sandro.bottaro@gmail.com',
      packages=["barnaba"],
      python_requires='>=3.6',
      install_requires=['numpy','scipy','mdtraj','future','scikit-learn'],
      scripts=['bin/barnaba'],
      zip_safe=False)


