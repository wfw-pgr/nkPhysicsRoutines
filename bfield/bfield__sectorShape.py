import os, sys, inspect
import numpy as np
import nkBasicAlgs.integrate__GaussLegendre as igl

# ============================================================ #
# =  integral formulae given by W.Kleeven ECPM (2016)        = #
# ============================================================ #

# ========================================================= #
# ===  bfield__sectorShape.py                           === #
# ========================================================= #

def bfield__sectorShape( coordinate=None, nGauss=15, bsign="+", zsign="+", \
                         r1    =None, r2    =None, th1=None, th2=None, z1=None, z2=None, \
                         alpha1=None, alpha2=None, J0 =None, \
                         coil_center=[0.,0.,0.], function__coilShape=None, coilkw={}, \
                         degree=True, append_coordinate=True ):

    x_, y_, z_ = 0, 1, 2

    # ------------------------------------------------- #
    # --- [1] arguments                             --- #
    # ------------------------------------------------- #
    if ( coordinate is None ):
        sys.exit( "[bfield__sectorShape.py]     coordinate == ??? [ERROR]" )
    if ( ( r1  is None ) or ( r2  is None ) ):
        sys.exit( "[bfield__sectorShape.py]     r1, r2     == ??? [ERROR]" )
    if ( ( th1 is None ) or ( th2 is None ) ):
        sys.exit( "[bfield__sectorShape.py]    th1, th2    == ??? [ERROR]" )
    if ( ( z1  is None ) or ( z2  is None ) ):
        sys.exit( "[bfield__sectorShape.py]     z1, z2     == ??? [ERROR]" )
    if ( ( alpha1 is None ) or ( alpha1 is None ) ):
        sys.exit( "[bfield__sectorShape.py] alpha1, alpha2 == ??? [ERROR]" )
    if ( J0 is None ):
        sys.exit( "[bfield__sectorShape.py]             J0 == ??? [ERROR]" )
    if ( degree ):
        th1   , th2    =    th1/180.0*np.pi,    th2/180.0*np.pi
        alpha1, alpha2 = alpha1/180.0*np.pi, alpha2/180.0*np.pi
    if ( zsign == "+" ):
        z1, z2 =  z1,  z2
    if ( zsign == "-" ):
        z1, z2 = -z1, -z2

    # ------------------------------------------------- #
    # --- [2] coil shape function :: z(r,th)        --- #
    # ------------------------------------------------- #
    def zCoil__linear( rp, th, zval=0.0, alpha=0.0 ):
        rM = 0.5 * ( np.min( rp ) + np.max( rp ) )
        return( zval + ( rp-rM ) * np.tan( alpha ) )
        
    if ( function__coilShape is None ):
        function__coilShape = zCoil__linear
    if ( not( inspect.isfunction( function__coilShape ) ) ):
        sys.exit( "[bfield__sectorShape.py] function__coilShape == ??? [ERROR]" )

    # ------------------------------------------------- #
    # --- [3] diffrential form :: function          --- #
    # ------------------------------------------------- #
    def Kleeven2016( rp, th, r0=0.0, th0=0.0, zval=0.0, alpha=0.0, J0=1.0, \
                     zCoil=zCoil__linear, coilkw={} ):
        zp     = zCoil( rp,th, zval=zval, alpha=alpha, **coilkw )
        insqrt = np.sqrt( + ( rp - r0*np.cos( th-th0 ) )**2 \
                          + (      r0*np.sin( th-th0 ) )**2 + zp**2 )
        bvalue = J0/( 4.0*np.pi ) * zp*rp * np.cos( alpha ) / insqrt**3
        return( bvalue )
        
    # ------------------------------------------------- #
    # --- [4] parameters for calculation            --- #
    # ------------------------------------------------- #
    x1Range  = [  r1,  r2 ]
    x2Range  = [ th1, th2 ]
    r0       = np.sqrt( + ( coordinate[:,x_]-coil_center[x_] )**2 \
                        + ( coordinate[:,y_]-coil_center[y_] )**2 )
    th0      = np.arctan2( coordinate[:,y_]-coil_center[y_], \
                           coordinate[:,x_]-coil_center[x_] )
    nPoints  = coordinate.shape[0]
    bfield   = np.zeros( (nPoints,3) )
    bz1,bz2  = np.zeros( (nPoints,) ), np.zeros( (nPoints,) )
    
    # ------------------------------------------------- #
    # --- [5] double integrate                      --- #
    # ------------------------------------------------- #
    for ik in list( range( nPoints ) ):
        kwargs1  = { "r0":r0[ik], "th0":th0[ik], "zCoil":function__coilShape, \
                     "J0":J0, "zval":z1, "alpha":alpha1, "coilkw":coilkw }
        kwargs2  = { "r0":r0[ik], "th0":th0[ik], "zCoil":function__coilShape, \
                     "J0":J0, "zval":z2, "alpha":alpha2, "coilkw":coilkw }
        bz1[ik]  = igl.integrate__GaussLegendre( nGauss =nGauss , function=Kleeven2016, \
                                                 x1Range=x1Range, x2Range=x2Range, \
                                                 kwargs =kwargs1 )
        bz2[ik]  = igl.integrate__GaussLegendre( nGauss =nGauss , function=Kleeven2016, \
                                                 x1Range=x1Range, x2Range=x2Range, \
                                                 kwargs =kwargs2 )
    if ( ( bsign == "+" ) and ( zsign == "+" ) ):
        bfield[:,z_] = + bz1[:] - bz2[:]
    if ( ( bsign == "-" ) and ( zsign == "+" ) ):
        bfield[:,z_] = - bz1[:] + bz2[:]
    if ( ( bsign == "+" ) and ( zsign == "-" ) ):
        bfield[:,z_] = - bz1[:] + bz2[:]
    if ( ( bsign == "-" ) and ( zsign == "-" ) ):
        bfield[:,z_] = + bz1[:] - bz2[:]
    
    # ------------------------------------------------- #
    # --- [6] return                                --- #
    # ------------------------------------------------- #
    if ( append_coordinate ):
        bfield = np.concatenate( [coordinate,bfield], axis=1 )
    return( bfield )


# ========================================================= #
# ===   Execution of Pragram                            === #
# ========================================================= #

if ( __name__=="__main__" ):

    x_, y_, z_ = 0, 1, 2
    bx_,by_,bz_= 3, 4, 5

    # ------------------------------------------------- #
    # --- [1] parameters                            --- #
    # ------------------------------------------------- #
    r1    , r2     =   0.4,  0.5
    th1   , th2    = -30.0, 30.0
    z1    , z2     =  0.02, 0.03
    alpha1, alpha2 =   0.0,  0.0
    J0             =   2.05
    nGauss         =   20
    
    # ------------------------------------------------- #
    # --- [2] generate coordinate                   --- #
    # ------------------------------------------------- #
    import nkUtilities.equiSpaceGrid as esg
    x1MinMaxNum = [   0.0,  1.0, 101 ]
    x2MinMaxNum = [  -0.5, +0.5, 101 ]
    x3MinMaxNum = [   0.0,  0.0,   1 ]
    coordinate  = esg.equiSpaceGrid( x1MinMaxNum=x1MinMaxNum, x2MinMaxNum=x2MinMaxNum, \
                                     x3MinMaxNum=x3MinMaxNum, returnType = "point"     )

    # ------------------------------------------------- #
    # --- [3] calculate bfield                      --- #
    # ------------------------------------------------- #
    bfield1 = bfield__sectorShape( coordinate=coordinate, nGauss=nGauss, bsign="+", zsign="+", \
                                   r1=r1, r2=r2, th1=th1, th2=th2, z1=z1, z2=z2, \
                                   alpha1=alpha1, alpha2=alpha2,   J0=J0, \
                                   append_coordinate=True )
    bfield2 = bfield__sectorShape( coordinate=coordinate, nGauss=nGauss, bsign="+", zsign="-", \
                                   r1=r1, r2=r2, th1=th1, th2=th2, z1=z1, z2=z2, \
                                   alpha1=alpha1, alpha2=alpha2,   J0=J0, \
                                   append_coordinate=True )
    bfield        = np.copy( bfield1 )
    bfield[:,bz_] = bfield1[:,bz_] + bfield2[:,bz_]
    
    # ------------------------------------------------- #
    # --- [4] save in a file                        --- #
    # ------------------------------------------------- #
    import nkUtilities.save__pointFile as spf
    outFile   = "test/bfield.dat"
    spf.save__pointFile( outFile=outFile, Data=bfield )

    # ------------------------------------------------- #
    # --- [5] draw colormap                         --- #
    # ------------------------------------------------- #
    import nkUtilities.load__config   as lcf
    import nkUtilities.cMapTri        as cmt
    config                   = lcf.load__config()
    pngFile                  = "test/bfield.png"
    config["FigSize"]        = [6,6]
    config["cmp_position"]   = [0.16,0.16,0.90,0.86]
    config["xTitle"]         = "X (m)"
    config["yTitle"]         = "Y (m)"
    config["cmp_xAutoRange"] = True
    config["cmp_yAutoRange"] = True
    config["cmp_xRange"]     = [-1.0,+1.0]
    config["cmp_yRange"]     = [-1.0,+1.0]

    cmt.cMapTri( xAxis=bfield[:,x_], yAxis=bfield[:,y_], cMap=bfield[:,bz_], \
    	         pngFile=pngFile, config=config )

    
