from setuptools import setup

setup(name='barnaba',
      version='0.1',
      description='analyze nucleic acid 3D structures and MD trajectories',
      url='https://github.com/srnas/barnaba',
      author='Sandro Bottaro',
      author_email='sandro dot bottaro at gmail dot com',
      packages=['barnaba'],
      install_requires=['python>=2.7','numpy','scipy'],
      test_suite='nose.collector',
      zip_safe=False)

