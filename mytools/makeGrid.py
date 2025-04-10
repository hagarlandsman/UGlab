
#I accidently deleted it...
#latitude,longitude
#33.089999999999996,35.449999999999996
#33.089999999999996,35.450399999999995


LatManara=33.1808 #33.194
LonManara=35.5458 #35.547
R=2 # in Km

LatToKm=110.574 # km/deg
LonToKm=111.320*0.837   # 111.3 * cos(Lat)

lat0=LatManara - R/LatToKm
lat1=LatManara + R/LatToKm
lon0=LonManara - R/LonToKm
lon1=LonManara + R/LonToKm
print (R,"km, is ",lat1-LatManara," Lon:",lon1-LonManara)
n=0
N=100
lonStep=(lon1-lon0)/N
latStep=(lat1-lat0)/N
print ("#N=",N,"  lon=",lon0,lon1,"\t lat=",lat0,lat1)
exit
lat=lat0


fileNumber=0
fName="ManaraGrid_"+str(R)+"km.N"+str(N)+"."+str(fileNumber)+".txt"
f=open(fName,"w")
f.write("latitude,longtitude\n")

while lat<=lat1:
    lon=lon0
    while lon<=lon1:
        if (n > 500*500):
            print ("n");
            n=0
            fileNumber=fileNumber+1
            fName="ManaraGrid_"+str(R)+"km.N"+str(N)+"."+str(fileNumber)+".txt"
            f.close()
            f=open(fName,"w")
            f.write("latitude,longtitude\n")

#        print (str(lat)+","+str(lon))
        f.write(str(lat)+","+str(lon)+"\n")
        lon=lon+lonStep
        n=n+1

    lat=lat+latStep
print ("n=",n)

