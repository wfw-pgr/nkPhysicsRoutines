import os, sys
import ctypes
import numpy   as np

# ================================================================ #
# ===  calc__biotsavartbfield                                  === #
# ================================================================ #

def calc__biotSavartBField( bfield=None, coils=None, I0=None ):
    
    # ---------------------------------------- #
    # --- [1]   引数チェック               --- #
    # ---------------------------------------- #
    if ( bfield is None ): sys.exit( "[calc__biotsavartbfield] bfield ???" )
    if ( coils  is None ): sys.exit( "[calc__biotsavartbfield] coils  ???" )
    if ( I0     is None ): sys.exit( "[calc__biotsavartbfield] I0     ???" )

    # ---------------------------------------- #
    # --- [2]   引数準備                   --- #
    # ---------------------------------------- #
    #  -- [2-1] 使用する引数を準備         --  #
    nBpt,nLpt = bfield.shape[0], coils.shape[0]
    #  -- [2-2] Fortranサイズへ変換        --  #
    bfield_   =     np.array( bfield, dtype=np.float64  )
    coils_    =     np.array( coils , dtype=np.float64  )
    I0_       =     np.array( I0    , dtype=np.float64  )
    nBpt_     = ctypes.byref( ctypes.c_int64( nBpt  )   )
    nLpt_     = ctypes.byref( ctypes.c_int64( nLpt  )   )

    # ---------------------------------------- #
    # --- [3]   ライブラリをロード         --- #
    # ---------------------------------------- #
    #  -- [3-1] ライブラリを定義           --  #
    path   = os.path.dirname( os.path.abspath( __file__ ) )
    pyLIB  = np.ctypeslib.load_library( 'pylib.so', path )
    #  -- [3-2] 入出力管理                 --  #
    pyLIB.calc__biotsavartbfield_.argtypes = [
        np.ctypeslib.ndpointer( dtype=np.float64 ),
        np.ctypeslib.ndpointer( dtype=np.float64 ),
        np.ctypeslib.ndpointer( dtype=np.float64 ),
        ctypes.POINTER( ctypes.c_int64   ),
        ctypes.POINTER( ctypes.c_int64   ),
    ]
    pyLIB.calc__biotsavartbfield_.restype = ctypes.c_void_p

    # ---------------------------------------- #
    # --- [4]   関数呼出 / 返却            --- #
    # ---------------------------------------- #
    pyLIB.calc__biotsavartbfield_( bfield_, coils_, I0_, nBpt_, nLpt_ )
    return( bfield_ )


# ================================================================ #
# ===  test execution                                          === #
# ================================================================ #
if ( __name__=='__main__' ):

    x_,y_,z_    = 0, 1, 2

    I0          = 100.0 * 1e3
    x0, y0, z0  = 0.0, 0.0, 1.0
    radius      = 1.0
    
    nLpt        = 101
    theta       = np.linspace( 0.0, 2.0*np.pi, nLpt )
    coils       = np.zeros( (nLpt,3) )

    coils[:,x_] = radius * np.cos( theta ) + x0
    coils[:,y_] = radius * np.sin( theta ) + y0
    coils[:,z_] =                          + z0
    
    import nkUtilities.equiSpaceGrid as esg
    x1MinMaxNum = [ -1.0, 1.0, 101 ]
    x2MinMaxNum = [ -1.0, 1.0, 101 ]
    x3MinMaxNum = [  0.0, 0.0,   1 ]
    ret         = esg.equiSpaceGrid( x1MinMaxNum=x1MinMaxNum, x2MinMaxNum=x2MinMaxNum, \
                                     x3MinMaxNum=x3MinMaxNum, returnType = "point" )
    nBpt        = ret.shape[0]
    bfield      = np.zeros( (nBpt,6) )
    bfield[:,x_:z_+1] = ret

    bfield      = calc__biotSavartBField( bfield, coils, I0 )
    print( bfield.shape )

    import nkUtilities.save__pointFile as spf
    outFile     = "dat/bfield_biot.dat"
    spf.save__pointFile( outFile=outFile, Data=bfield )
