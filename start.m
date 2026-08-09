function start()
%START Initialize the MATLAB ISAC library path for this session.

projectRoot = fileparts(mfilename("fullpath"));
sourceFolder = fullfile(projectRoot, "src");

if ~isfolder(sourceFolder)
    error("isac:startup:MissingSource", ...
        "ISAC source folder was not found: %s", sourceFolder);
end

addpath(sourceFolder);
fprintf("ISAC library initialized: %s\n", sourceFolder);
end
