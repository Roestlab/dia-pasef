from setuptools import setup

def readme():
    with open('README.rst') as f:
        return f.read()

setup(name='diapysef',
      version='0.3.3',
      description='Analysis and visualization of DIA PASEF data',
      long_description=readme(),
      author='Max Frank, Annie Ha, Hannes Roest',
      author_email='hannes.rost@utoronto.ca',
      license='MIT',
      packages=['diapysef'],
      install_requires=[
          'pandas',
          'numpy',
          'matplotlib',
          'statsmodels',
          'pyopenms',
          'patsy'],
      package_data={
          'diapysef': ['data/*']
      },
      scripts=['scripts/get_dia_windows.py',
               'scripts/annotate_mq_ionmobility.py',
               'scripts/plot_dia_windows.py',
               'scripts/convertTDFtoMzML.py'],
      zip_safe=False)
