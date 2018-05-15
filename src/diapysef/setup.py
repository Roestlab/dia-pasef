from setuptools import setup

def readme():
    with open('README.rst') as f:
        return f.read()

setup(name='diapysef',
      version='0.1',
      description='Analysis and visualization of DIA PASEF data',
      long_description=readme(),
      author='Max Frank, Hannes Roest',
      author_email='max.frank@mail.utoronto.ca',
      license='MIT',
      packages=['diapysef'],
      install_requires=[
          'pandas',
          'numpy',
          'matplotlib',
          'statsmodels'
          'pyopenms'],
      package_data={
          'diapysef': ['data/*']
      },
      scripts=['scripts/get_dia_windows.py',
               'scripts/annotate_mq_ionmobility.py',
               'scripts/plot_dia_windows.py',
               'scripts/convertTDFtoMzML.py'],
      zip_safe=False)
