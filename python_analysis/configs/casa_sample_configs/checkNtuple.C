bool checkNtuple(string ntuple, int n_core) {
	ROOT::EnableImplicitMT(n_core);
	ROOT::TThreadedObject <TH1F> Bad("Bad", "Bad", 1,0,1);
	ROOT::TTreeProcessorMT tp(ntuple.c_str(), "ntuples/outT");
	auto myFunction = [&](TTreeReader &myReader) {
		auto bad = Bad.Get();
		while (myReader.Next()) {
			if (myReader.GetEntryStatus() > 0) {
				cout << "Bad entry status " << myReader.GetEntryStatus() << endl;
				bad->Fill(0.5);
				break;
			}
		}
		if (myReader.GetEntryStatus() != TTreeReader::kEntryBeyondEnd) {
			cout << "Filling at end with " << myReader.GetEntryStatus() << " vs " << TTreeReader::kEntryBeyondEnd << endl;
			bad->Fill(0.5);
		}
   };
   tp.Process(myFunction);
   auto mergeBad = Bad.Merge();
   cout << mergeBad->Integral() << endl;
   if (mergeBad->Integral() > 0) {
		return true;
   }
   else {
		return false;
   }
}