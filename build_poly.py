import numpy as np
import h5py

def makePolygons(mask):
    """Return a set of polygons from the mask bounding the >0 pixels

    Returns:
    [ [(x11,y11), (x12,y12)...], [[(x21,y21),...]] ]
    """

    # unprocessed pixels are set to be -1, so we can keep track of the
    # output polygons as >1
    msk = np.where(mask>0, -1, 0)

    # Make a list of line segments around a pixel
    # Each segment consists of (initial_y, initial_x, direction)

    # DIRECTIONS
    #   1 ->  
    #   ----- 2
    #/\ |   | |
    # | |   | V 
    # 4 |   |
    #   -----
    #   <- 3

    # If a segment lies next to an unprocessed mask pixel, we can
    # replace the segment with three which also encompass the pixel

    RIGHT, DOWN, LEFT, UP = 1, 2, 3, 4
    deltadirn = {RIGHT: (0,1), DOWN: (-1,0), LEFT: (0,-1), UP: (1,0)}
    yw, xw = msk.shape

    polys = []
    polyidx = 1   # keep track of which pixels assigned to which poly
    while True:
        # find the next unprocessed pixel to start
        idx0 = np.where(msk<0)
        if len(idx0[0]) == 0:
            # all done, so exit
            break
        y0, x0 = idx0[0][0], idx0[1][0]

        # initial segments around starting pixel
        segs = [
            (y0, x0, UP), (y0+1,x0, RIGHT),
            (y0+1,x0+1, DOWN), (y0,x0+1, LEFT)]
        msk[y0,x0] = polyidx

        iseg = 0
        while iseg < len(segs):
            y0, x0, dirn = segs[iseg]

            newseg = None
            if dirn == RIGHT:
                # segment is rightward (+x)
                if y0 < yw and msk[y0,x0]<0:
                    newseg = [(y0,x0,UP), (y0+1,x0,RIGHT), (y0+1,x0+1,DOWN)]
                    msk[y0,x0] = polyidx
                elif y0 > 0 and msk[y0-1,x0]<0:
                    newseg = [(y0,x0,DOWN), (y0-1,x0,RIGHT), (y0-1,x0+1,UP)]
                    msk[y0-1,x0] = polyidx
            elif dirn == DOWN:
                # segment is downward (-y)
                if x0 < xw and msk[y0-1,x0]<0:
                    newseg = [(y0,x0,RIGHT), (y0,x0+1,DOWN), (y0-1,x0+1,LEFT)]
                    msk[y0-1,x0] = polyidx
                elif x0 > 0 and msk[y0-1,x0-1]<0:
                    newseg = [(y0,x0,LEFT), (y0,x0-1,DOWN), (y0-1,x0-1,RIGHT)]
                    msk[y0-1,x0-1] = polyidx
            elif dirn == LEFT:
                # segment is leftward (-x)
                if y0 < yw and msk[y0,x0-1]<0:
                    newseg = [(y0,x0,UP), (y0+1,x0,LEFT), (y0+1,x0-1,DOWN)]
                    msk[y0,x0-1] = polyidx
                elif y0 > 0 and msk[y0-1,x0-1]<0:
                    newseg = [(y0,x0,DOWN), (y0-1,x0,LEFT), (y0-1,x0-1,UP)]
                    msk[y0-1,x0-1] = polyidx
            elif dirn == UP:
                # segment is upward (+y)
                if x0 > 0 and msk[y0,x0-1]<0:
                    newseg = [(y0,x0,LEFT), (y0,x0-1,UP), (y0+1,x0-1,RIGHT)]
                    msk[y0,x0-1] = polyidx
                elif x0 < xw and msk[y0,x0]<0:
                    newseg = [(y0,x0,RIGHT), (y0,x0+1,UP), (y0+1,x0+1,LEFT)]
                    msk[y0,x0] = polyidx

            if newseg is not None:
                # replace segment with new segments
                segs = segs[:iseg] + newseg + segs[iseg+1:]
            else:
                # move to next otherwise
                iseg += 1

        # clean up line segments where we go in one direction, then
        # back again
        i = 0
        while i<len(segs):
            d1 = segs[i][2]
            d2 = segs[(i+1)%len(segs)][2]
            if( (d1==LEFT and d2==RIGHT) or (d1==RIGHT and d2==LEFT) or
                (d1==UP and d2==DOWN) or (d1==DOWN and d2==UP) ):
                del segs[i]
                del segs[i % len(segs)]
                if i > 0:
                    i -= 1
            else:
                i += 1

        # move any with the same direction at end to start to
        # maximimize runs
        while segs[0][2] == segs[-1][2]:
            segs.insert(0, segs.pop())

        # turn line segments into a polygon
        poly = []
        lastdir = None
        for seg in segs:
            d = deltadirn[seg[2]]
            endpt = (seg[0]+d[0], seg[1]+d[1])
            if lastdir == seg[2]:
                # if we go in the same direction, replace last point
                poly[-1] = endpt
            else:
                poly.append(endpt)
            lastdir = seg[2]

        polys.append(poly)
        polyidx += 1

    print(msk)
    return polys

def main():
    mask = np.array([
        [0,1,1,0,0,1],
        [1,1,1,1,0,1],
        [0,1,0,1,0,0],
        [0,1,1,1,0,0],
        [0,0,0,0,1,0],
    ], dtype=np.int16)

    polys = makePolygons(mask)

    with h5py.File('out.hdf5','w') as fout:
        for i, poly in enumerate(polys):
            for p in poly:
                print('%g,%g' % tuple(p))
            print()

            poly.insert(0, poly[-1])
            fout[f'x{i+1}'] = np.column_stack(poly)[0,:]
            fout[f'y{i+1}'] = np.column_stack(poly)[1,:]

if __name__ == '__main__':
    main()
