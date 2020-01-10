from setuptools import setup

def readme():
    with open('README.rst') as f:
        return f.read()

desc = """\
Analysis, conversion and visualization of diaPASEF data."""

setup(name='diapysef',
      version='0.3.5',
      description=desc,
      long_description=readme(),
      url="https://github.com/Roestlab/dia-pasef",
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
               'scripts/convertTDFtoMzML.py',
               'scripts/create_library.py',
               'scripts/high_precision_irt.py'],
      zip_safe=False)
