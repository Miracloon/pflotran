rm -rf output/
python3 create_dfn.py
python3 ./mapdfn2pflotran.py ./output
