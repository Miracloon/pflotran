def get_dimensions():

    # for 543 irregular grid spacing, dx,dy,dz must be specified.
    #nx = 5; ny = 4; nz = 3
    #dx = [10.,11.,12.,13.,14.]
    #dy = [13.,12.,11.,10.]
    #dz = [15.,20.,25.]

    # for uniform grid spacing, specify one value in each direction
    nx = 40; ny = 40; nz = 20
    dx = [50.]
    dy = [50.]
    dz = [10.]

    # note that you can mix and match irregular spacing along the different axes
    n = (nx,ny,nz)
    d = [dx,dy,dz]
    origin = (0.,0.,0.)

    return n, d, origin
