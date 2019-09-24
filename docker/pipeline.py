import papermill as pm

pm.execute_notebook(
   'Pipeline.ipynb',
   'Pipeline_out.ipynb',
   parameters = dict(
   workingdir="/home/filippo/tacos/files/",
   execdir="/home/filippo/tacos/")
)
