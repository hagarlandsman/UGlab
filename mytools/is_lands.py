from global_land_mask import globe

fname = "R100km_0.1.csv"
fname = "R10km_0.1.csv"
with open(fname+".out","w") as fo:
    with open(fname) as f:
        for l in f:
            isl = -1
           # print (l)
            if (l.find("lat")>0):
                next
            a=l.split(",")

            #print (a[0],a[1])
            try:
                if globe.is_land(float(a[0]),float(a[1])):
                    isl = 1
                else:
                    isl = 0
            except Exception  as e:
                print (e)
            fo.write(l.strip()+","+str(isl)+"\n")