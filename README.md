Traceback (most recent call last):
  File "/insaflu_web/TELEVIR/main.py", line 215, in <module>
    main()
  File "/insaflu_web/TELEVIR/main.py", line 180, in main
    metagen_prep.setup_soft()
  File "/insaflu_web/TELEVIR/install_scripts/main_install.py", line 439, in setup_soft
    self.dl_metadata_prot()
  File "/insaflu_web/TELEVIR/install_scripts/main_install.py", line 295, in dl_metadata_prot
    self.wdir.generate_main_protacc_to_taxid()
  File "/insaflu_web/TELEVIR/install_scripts/modules/db_install.py", line 490, in generate_main_protacc_to_taxid
    general_db = pd.concat(to_concat, axis=0)
  File "/opt/venv/lib/python3.9/site-packages/pandas/util/_decorators.py", line 311, in wrapper
    return func(*args, **kwargs)
  File "/opt/venv/lib/python3.9/site-packages/pandas/core/reshape/concat.py", line 347, in concat
    op = _Concatenator(
  File "/opt/venv/lib/python3.9/site-packages/pandas/core/reshape/concat.py", line 404, in __init__
    raise ValueError("No objects to concatenate")
ValueError: No objects to concatenate
ERROR: Service 'televir-server' failed to build : The command '/bin/sh -c /opt/venv/bin/python main.py --docker --seqdl --partial' returned a non-zero code: 1