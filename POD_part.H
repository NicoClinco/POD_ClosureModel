// Creating the PtrList for fields and modes

PtrList<volVectorField> Vfield; 
PtrList<volScalarField> Sfield;

PtrList<volVectorField> Vmodes;
PtrList<volScalarField> Smodes;


runTime.setTime(Times[startTime], startTime);
 
wordList fieldlist
(
 ITHACAPODdict.lookup("fields")
);
 
 
if (startTime == endTime)
{
    Info << "The case has no snapshots to process, exiting the code" << endl;
    exit(0);
}

//  For all the fields ( eg. U and p ) do a loop 

for (label k = 0; k < fieldlist.size(); k++)
{
       dictionary& subDict = ITHACAPODdict.subDict(fieldlist[k]);
       scalar nmodes = readScalar(subDict.lookup("nmodes"));
       word field_name(subDict.lookup("field_name"));
       word field_type(subDict.lookup("field_type"));

 //  Perform the real POD // 
	for (label i = startTime; i < endTime + 1; i++)
	  {
            Info << "Reading snapshot " << i << " for field " << field_name << endl;
            runTime.setTime(Times[i], i);
            mesh.readUpdate();
	    if (field_type == "vector")
	    {
            volVectorField vector_field
            (
                    IOobject
                    (
                        field_name,
                        runTime.timeName(),
                        mesh,
                        IOobject::MUST_READ
                    ),
                    mesh
             );
	    Vfield.append(vector_field.clone());

	    }else{
                volScalarField scalar_field
                (
                    IOobject
                    (
                        field_name,
                        runTime.timeName(),
                        mesh,
                        IOobject::MUST_READ
                    ),
                    mesh
                );
                Sfield.append(scalar_field.clone());

	    }// end switch scalar vs vector

	  }// End fields-acquisition

   /* ------- Get the modes ------------ */ 	
    
	    if (field_type =="vector")
	    {
		ITHACAPOD::getModes(Vfield, Vmodes, field_name, 0, 0, 0, nmodes);
            }
	    if (field_type == "scalar")
	    {
		 ITHACAPOD::getModes(Sfield, Smodes, field_name, 0, 0, 0, nmodes);
	    }
  /* ------------------------------------ */
}//  end Function
