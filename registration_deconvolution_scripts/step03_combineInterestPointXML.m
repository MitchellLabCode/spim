% Script to combine xml files with different interest points recorded
% After detecting interest points from different channels/timepoints in
% parallel across different Fiji instances, this code combines these
% detections into a single XML file.

xmlfns = {'./dataset_ch0IPs.xml', ...
    './dataset_ch1IPs.xml', ...
    './dataset_ch2IPs.xml'} ;
outfn = './dataset_allIPs.xml' ;
combine_xml_files(xmlfns, outfn)