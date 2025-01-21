# We need a python function that when called would produce a gif
# - grab the profile
# - grab the mmcif
# - open mmcif in chimerax
# - paint it with ribrep
# - record a movie
# - on completiong: ffmpeg it to a gif


RCSB_ID = "7K00"
def produce_gif(rcsb_id:str):
    render_movie = 'chimerax --nogui --offscreen --cmd "cd /home/rtviii/dev/riboxyz/ribctl/lib/chimerax; open cmd_ribetl.py; open cmd_ribrepr.py; open cmd_ribmovie.py; ribmovie {}; close all" -'.format(rcsb_id, rcsb_id)
    import subprocess
    subprocess.run(render_movie, shell=True)

produce_gif(RCSB_ID)

