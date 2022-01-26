function figSave(filename,filetype,figFlag)

if figFlag == 1
    export_fig([filename,filetype]);
    savefig([filename,'.fig']);
else
        export_fig([filename,filetype]);
end
close all force;