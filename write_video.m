function write_video(directory,video_name)
PNG_all = fullfile(directory, '*.PNG');
PNG_dir = dir(PNG_all);
video_obj = VideoWriter(video_name);
video_obj.FrameRate = 3;
open(video_obj);

for num_frame = 1:length(PNG_dir)
    image_name = PNG_dir(num_frame).name;
    image_dir = fullfile(directory,image_name);
%     display image name in the command window
    fprintf(1, 'Now reading %s\n', image_name);
%     Write this frame out to the AVI file.
    writeVideo(video_obj, imread(image_dir));
end
close(video_obj);
end

