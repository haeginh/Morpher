#ifndef README_HH
#define README_HH
#include <iostream>

void PrintUsage(){
    cout<<"<Morphing Procedure>"<<endl<<endl;

    cout<<"-preprocessing-"<<endl;
    cout<<"1. Prepare OBJ and MESH file for whole phantom."<<endl;
    cout<<"2. Prepare precise frame phantom (MESH) to fully cover the whole phantom."<<endl;
    cout<<"3. Calculate BBW for bones and joints."<<endl;
    cout<<"\t./Morpher -bbw [PREFIX.mesh] [frame.ply] [tgf file] [bone P interval(default:1cm)]"<<endl;
    cout<<"\toutput: PREFIX.fbary, PREFIX_bbwF.mesh, PREFIX_fb.w, PREFIX_fs.w, PREFIX.tgf"<<endl<<endl;
    cout<<"4. Widen the arms and legs of the whole phantom."<<endl;
    cout<<"\t./Morpher -bary [PREFIX(.mesh)] [shell.obj]"<<endl;
    cout<<"\toutput: PREFIX.bary, PREFIX_s.mesh, PREFIX_s.fbary"<<endl<<endl;

    cout<<"-morphing-"<<endl;
    cout<<"5. Create frame phantom in RapidForm (decimate(vertice preserve)) and save as filename_frame.ply"<<endl;
    cout<<"6. Calculate barycentric coord."<<endl;
    cout<<"\t./Morpher [PREFIX] [target.obj]"<<endl;
    cout<<"\toutput: filename_frame.mesh, filename_bary.m"<<endl;

    cout<<"-morphing-"<<endl;
}

#endif // README_HH
