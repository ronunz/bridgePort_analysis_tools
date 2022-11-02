void bridgeport_spectrum_histo_generator()
{
    TH1F *histogram = new TH1F("Histogram",1000,0,1024);

    fstream file;
    file.open("Cs-137_1_m_5_min.csv", ios::in);

    double numberOfCountsInChannel;

    while(i)
    {
        file >> numberOfCountsInChannel;
        int channelNumber;
        for (int i=0; i < numberOfCountsInChannel; i++)
        {
            histogram->Fill(channelNumber);
        }
    }

    file.close();
}