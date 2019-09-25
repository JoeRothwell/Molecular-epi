# Instructions for running circos on PC
# Install perl 5.16 and circos into C:
# Check perl is installed
perl -v
# Try and run circos (command prompt not git bash)
perl C:/circos/bin/circos -version
# Navigate to Perl directory and run cpan
cpan
# Install missing modules (look at list)

install Pod::Usage
install Font::TTF:Font
install Config::General
install SVG
install Math::Bezier
install Params:Validate
install Math::Round
install Math::VecStat
install Regexp::Common
install Statistics::Basic
install Text::Format
install Set::IntSpan
exit

# Put karyotype, data files and .conf file in C:/circos and run
perl C:/circos/bin/circos -conf test.conf

# output is found in C:/perl/bin