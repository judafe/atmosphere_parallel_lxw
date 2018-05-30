clear
% Read data from .txt files or load previous run data. ----------
readfile01 = 1;
% 1 -- > read from .txt files
% else > load previous run data
% ----------------------------------------------------------------
loadrun = 'run2';
prefix = {' rho',' rhoU',' rhoV',' energy'};


if readfile01 == 1 
    xmin = -100000; xmax = 100000;
    ymin = 0; ymax = 160000;
    dx = 200; dy = 200;
    
    xcv=xmin-3*dx/2:dx:xmax+3*dx/2;
    ycv=ymin-3*dy/2:dy:ymax+3*dy/2;
    
    xdomain = xmax-xmin;
    ydomain = ymax-ymin;
    
    ncx = xdomain/dx+4; % number of cell centers (I in matlab code)
    ncy = ydomain/dy+4; % number of cell centers (J in matlab code)
    
    nprocs =4;
    ncxs = (ncx-4)/nprocs + 4;
    
    nsteps = 500;
    
    Qfiles = {-1*ones(1,nprocs)...
        -1*ones(1,nprocs)...
        -1*ones(1,nprocs)...
        -1*ones(1,nprocs)};
    
   
    for npter = 1:nprocs
        for qter = 1:4
            Qfiles{qter}(npter) = fopen([prefix{qter},num2str(npter-1),'.txt']);
        end
    end
    
    Qhist = {NaN*ones((ncy-4),(ncxs-4),nsteps)...
        NaN*ones((ncy-4),(ncxs-4),nsteps)...
        NaN*ones((ncy-4),(ncxs-4),nsteps)...
        NaN*ones((ncy-4),(ncxs-4),nsteps)};
    
    Qcol={[],[],[],[]};
    for kter = 1:4 % loop through Q
        Qcol{kter} = NaN*ones((ncy-4)*(ncxs-4)*nsteps,1);
    end
    
    
    for npter = 1:nprocs % loop through processors
        for qter = 1:4 % loop through Q
            % Qfiles{qter}(npter) = fopen([prefix{qter},num2str(npter-1),'.txt']);
            Qcol{qter} = fscanf(Qfiles{qter}(npter),'%f',(ncy-4)*(ncxs-4)*nsteps);
            colsize=size(Qcol{qter});
            disp(['Q(',num2str(qter),'), from file: ',[prefix{qter},num2str(npter-1),'.txt'],' ',...
                num2str(colsize)])
            if Qfiles{qter}(npter)==-1
                error([prefix{qter},num2str(npter-1),'.txt file not found'])
            elseif colsize(1)==0
                error([prefix{qter},num2str(npter-1),'.txt file data not loaded'])
            end
        end
        
        xstart = (ncxs-4)*(npter-1)+1;
        xfin = (ncxs-4)*npter;
        
        
        for iter = 1:nsteps % loop through recordings
            ystart = (ncy-4)*(iter-1)+1;
            yfin = (ncy-4)*iter;
            colstart = (iter-1)*(ncy-4)*(ncxs-4)+1;
            colfin   = iter*(ncy-4)*(ncxs-4);
            
            for qter = 1:4 % loop through Qs
                Qhist{qter}(:,xstart:xfin,iter) ...
                    = reshape(Qcol{qter}(colstart:colfin),[ncxs-4,ncy-4])';
            end
        end
    end
else
    load(loadrun)
end % end of read-file section

figure(1); clf
for iter = 1:nsteps
    for kter = 1:4
        subplot(2,2,kter)
        imagesc(xcv(3:end-2)/1000,ycv(3:end-2)/1000,...
            Qhist{kter}(:,:,iter)); colorbar
        axis xy
        title(prefix{kter})
        pause(0.02)
    end
end