pro gasiords_quies

; tracks which gas particles are each by which black holes, and how
; much is eaten from each.
; the goal is to trace the gas particles and find out whether the
; accretion is merger-induced or cold flow or what.

; output is a file for each BH, with the iords of each gas particle
; that hole accretes and the amount of mass (in Msun) accreted from
; each.


;readcol,'out.distance',bhsink,bhiord,w,time,w,w,w,w,w,gasiord,d,w,w,l,format='a,l,a,f,a,a,a,a,a,l,d,a,a,d',/silent

readcol,'out.dmq',bhsink,bhiord,gasiord,w,dmq,qcurlv,pcurlv,format='a,l,l,a,d,d,d',/silent
munit=1.84793e16

unique_iords = bhiord[uniq(bhiord,sort(bhiord))]

n=n_elements(unique_iords)
; number of bhs
print,n,' black holes'
close,1

for i=0,n-1 do begin
    thisbh = where(bhiord eq unique_iords[i])
    eatengas = gasiord[thisbh]
    eatenmass = dmq[thisbh]
    ; only unique iords here.
    uniquegas=uniq(eatengas,sort(eatengas)) ; indices
    eacheatengas=eatengas[uniquegas]
    filename='gasiords.'+trim(unique_iords[i])
    openw,1,filename
    for j=0L,n_elements(eacheatengas)-1 do begin
        eacheatenmass=eatenmass[where(eatengas eq eacheatengas[j])]
        printf,1,eacheatengas[j],total(eacheatenmass)*munit/2.0
        ; divide by two since there's a duplicate output of each.
    endfor
    close,1
endfor


spawn,'ls gasiords.* > gasiords.list'


end

