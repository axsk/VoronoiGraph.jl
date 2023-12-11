using CoverageTools;
nlines(dir) = [sum(x -> (!isnothing(x) && x .> 0), c.coverage) for c in CoverageTools.process_folder(dir)];