double precision     ::  temperature(4) !temperature, K
double precision     ::  pressure(6) !pressure, bar
character(LEN=8)     ::  strPress(6),strTemp(4) 
character(LEN=8)     ::  species(4) !species considered
double precision     ::  moleweight(4) !molecule weight 
integer		     ::  nTemp,nPress,nSpecies
COMMON /PTSpeciesMw/  temperature,pressure,species,moleweight, &
                      nTemp,nPress,nSpecies,strPress,strTemp
