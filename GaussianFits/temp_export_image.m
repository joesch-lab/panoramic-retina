%adhoc script for eexporting the image used for gprops fitting (they were
%not saved initially)s


files = dir("/nfs/scistore08/joeschgrp/dgupta/TER2/V1/combined_*.mat");

%%
count = 1;
for retina = files(2:6)'
    count = count +1;
    load(fullfile(retina.folder, retina.name));
    
    disp(fullfile(retina.folder, retina.name));
    
    output(count).gprops = RF_temp_export_image(combined_datasets.rf(1));
    
    clear combined_datasets;
end

%%

cnt = 1;
ot = 1;

for ret = 1:6
    gprops = output(ret).gprops;
    for i = 1:numel(gprops)
        if sel_rf(cnt)
            outg(ot) = gprops(i);
            ot = ot + 1;
        end
        cnt = cnt + 1;
    end
end

%%
figure;
for i = 1:100:10000
    imagesc(outg(i).RFaroundpeak, [-1,1]);
    colormap(redblue);
    input("");
end
