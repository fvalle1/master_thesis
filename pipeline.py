import papermill as pm

pm.execute_notebook(
   'Pipeline.ipynb',
   'Pipeline_out.ipynb',
   parameters = dict(workingdir=r"/Volumes/GoogleDrive/My Drive/tesi_magistrale/tesi/gtex/hsbm/miRNA")
)
