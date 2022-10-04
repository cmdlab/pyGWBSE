from setuptools import setup


setup(
   name='pyGWBSE',
   version='1.0',
   description='python workflow for GW-BSE calculation',
   author='Tathagata Biswas',
   author_email='tbiswas3@asu.edu',
   packages=['pyGWBSE'],
   install_requires=['FireWorks>=2.0.3', 'pymatgen>=2022.9.21',            
                          'atomate>=1.0.3'],   
   package_data={"pyGWBSE": ["inputset.yaml"]}

)
