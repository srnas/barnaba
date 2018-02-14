from setuptools import setup,find_packages


setup(name='barnaba',
      version='0.1',
      description='analyze nucleic acid 3D structures and MD trajectories',
      url='https://github.com/srnas/barnaba',
      author='Sandro Bottaro',
      author_email='sandro dot bottaro at gmail dot com',
      packages=find_packages(),
      python_requires='>=2.6',
      install_requires=['numpy','scipy','mdtraj','future'],
      test_suite='nose.collector',
      zip_safe=False)


