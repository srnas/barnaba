from setuptools import setup,find_packages

def readme():
    with open('README.rst') as f:
        return f.read()
    
setup(name='barnaba',
      description='analyze nucleic acid 3D structures and MD trajectories',
      long_description=readme(),
      classifiers=[
          'Development Status :: 5 - Production/Stable',
          'Intended Audience :: Science/Research',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
          'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
          'Programming Language :: Python :: 2',
          'Programming Language :: Python :: 2.6',
          'Programming Language :: Python :: 2.7',
          'Programming Language :: Python :: 3',
          'Programming Language :: Python :: 3.3',
          'Programming Language :: Python :: 3.4',
          'Programming Language :: Python :: 3.5',
          'Programming Language :: Python :: 3.6',
      ],      
      url='https://github.com/srnas/barnaba',
      author='Sandro Bottaro',
      author_email='sandro.bottaro@gmail.com',
      use_scm_version = True,
      setup_requires = ['setuptools_scm'],
      packages=find_packages(),
      python_requires='>=2.6',
      install_requires=['numpy','scipy','mdtraj','future'],
      test_suite='nose.collector',
      scripts=['bin/barnaba'],
      zip_safe=False)


