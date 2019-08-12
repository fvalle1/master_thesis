import papermill as pm

pm.execute_notebook(
   'Pipeline.ipynb',
   'Pipeline_out.ipynb',
   parameters = dict(alpha=0.6, ratio=0.1)
)
