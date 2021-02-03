import pathlib
from setuptools import setup, find_packages

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

HERE = pathlib.Path(__file__).parent
README = (HERE / "README.md").read_text()

setup(name='scma',
      version='0.1.1',
      packages=find_packages(),
      description='Single-cell metabolomics analysis',
      long_description=README,

      author='Kang Tian',
      author_email='tiank18@mails.tsinghua.edu.cn',
      url='https://github.com/tkisss/metabolomics',
      scripts=['scMA.py'],
      install_requires=requirements,
      python_requires='>3.6.0',
      license='MIT',

      classifiers=[
          'Development Status :: 4 - Beta',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: MIT License',
          'Programming Language :: Python :: 3.7',
          'Operating System :: MacOS :: MacOS X',
          'Operating System :: Microsoft :: Windows',
          'Operating System :: POSIX :: Linux',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
     ],
     )
