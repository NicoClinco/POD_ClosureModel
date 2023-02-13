    const entry* existnsnap = ITHACAPODdict.lookupEntryPtr("Nsnapshots", false, true);
    const entry* existLT = ITHACAPODdict.lookupEntryPtr("FinalTime", false, true);
 
    // Initiate variable from PODSolverDict
    if ((existnsnap) && (existLT))
    {
        Info << "Error you cannot define LatestTime and NSnapShots together" << endl;
        abort();
    }
    else if (existnsnap)
    {
        scalar InitialTime = ITHACAPODdict.lookupOrDefault<scalar>("InitialTime", 0);
        nSnapshots = readScalar(ITHACAPODdict.lookup("Nsnapshots"));
        startTime = Time::findClosestTimeIndex(runTime.times(), InitialTime);
        nSnapshots = min(nSnapshots , Times.size() - startTime);
        endTime = startTime + nSnapshots - 1;
        Info <<"Selected : "<<  nSnapshots << "snapshots"<<endl;
    }
    else
    {
        scalar InitialTime = ITHACAPODdict.lookupOrDefault<scalar>("InitialTime", 0);
       // scalar OFfinalTime = runTime.controlDict().lookup("endTime");
        scalar FinalTime = ITHACAPODdict.lookupOrDefault<scalar>("FinalTime",1e6);
        endTime = Time::findClosestTimeIndex(runTime.times(), FinalTime);
        startTime = Time::findClosestTimeIndex(runTime.times(), InitialTime);
        nSnapshots = endTime - startTime + 1;
        if (InitialTime > FinalTime)
        {
            Info << "FinalTime cannot be smaller than the InitialTime check your ITHACAPODdict file\n" << endl;
            abort();
        }
    }
    // Print out some Infos
    Info << "startTime: " << startTime << "\n" << "endTime: " << endTime << "\n" << "nSnapshots: " << nSnapshots << "\n" << endl;
