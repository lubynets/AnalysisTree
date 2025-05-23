//
// Created by oleksii on 09.04.25.
//

#include "GenericContainerFiller.hpp"

using namespace AnalysisTree;

GenericContainerFiller::GenericContainerFiller(std::string fileInName, std::string treeInName) : file_in_name_(std::move(fileInName)),
                                                                                                 tree_in_name_(std::move(treeInName)) {}

void GenericContainerFiller::Init() {
  file_in_ = TFile::Open(file_in_name_.c_str(), "read");
  if (file_in_ == nullptr) throw std::runtime_error("GenericContainerFiller::Run(): file_in_ == nullptr");

  tree_in_ = file_in_->Get<TTree>(tree_in_name_.c_str());
  if (tree_in_ == nullptr) throw std::runtime_error("GenericContainerFiller::Run(): tree_in_ == nullptr");

  if (!fields_to_ignore_.empty() && !fields_to_preserve_.empty()) throw std::runtime_error("GenericContainerFiller::Run(): !fields_to_ignore_.empty() && !fields_to_preserve_.empty()");

  BranchConfig branchConfig(branch_out_name_, DetType::kGeneric);

  auto lol = tree_in_->GetListOfLeaves();
  const int nLeaves = lol->GetEntries();
  for (int iLeave = 0; iLeave < nLeaves; iLeave++) {
    auto leave = lol->At(iLeave);
    const std::string fieldName = leave->GetName();
    const std::string fieldType = leave->ClassName();
    if (!fields_to_ignore_.empty() && (std::find(fields_to_ignore_.begin(), fields_to_ignore_.end(), fieldName) != fields_to_ignore_.end())) continue;
    if (!fields_to_preserve_.empty() && (std::find(fields_to_preserve_.begin(), fields_to_preserve_.end(), fieldName) == fields_to_preserve_.end())) continue;
    if (fieldType == "TLeafF") {
      branchConfig.AddField<float>(fieldName);
    } else if (fieldType == "TLeafI" || fieldType == "TLeafB" || fieldType == "TLeafS") {
      branchConfig.AddField<int>(fieldName);
    }
    branch_map_.emplace_back((IndexMap){fieldName, fieldType, branchConfig.GetFieldId(fieldName)});
  }
  branch_values_.resize(branch_map_.size());

  config_.AddBranchConfig(branchConfig);

  for (int iV = 0; iV < branch_values_.size(); iV++) {
    TBranch* branch = tree_in_->GetBranch(branch_map_.at(iV).name_.c_str());
    SetAddressFICS(branch, branch_map_.at(iV), branch_values_.at(iV));
  }

  generic_detector_ = new GenericDetector(branchConfig.GetId());

  file_out_ = TFile::Open(file_out_name_.c_str(), "recreate");
  tree_out_ = new TTree(tree_out_name_.c_str(), "Analysis Tree");
  tree_out_->SetAutoSave(0);
  tree_out_->Branch((branchConfig.GetName() + ".").c_str(), "AnalysisTree::GenericDetector", &generic_detector_);

  entry_switch_trigger_id_ = entry_switch_trigger_var_name_.empty() ? -999 : DetermineFieldIdByName(branch_map_, entry_switch_trigger_var_name_);
}

int GenericContainerFiller::Exec(int iEntry, int previousTriggerVar) {
  tree_in_->GetEntry(iEntry);
  const int currentTriggerVar = entry_switch_trigger_id_ >= 0 ? branch_values_.at(entry_switch_trigger_id_).get() : previousTriggerVar;
  auto isNewATEntry = [&]() { return iEntry == 0 || (currentTriggerVar != previousTriggerVar) || (n_channels_per_entry_ >= 0 && iEntry % n_channels_per_entry_ == 0); };

  if (isNewATEntry()) {
    if (iEntry != 0) tree_out_->Fill();
    generic_detector_->ClearChannels();
  }
  auto& channel = generic_detector_->AddChannel(config_.GetBranchConfig(generic_detector_->GetId()));
  SetFieldsFICS(branch_map_, channel, branch_values_);

  return currentTriggerVar;
}

void GenericContainerFiller::Finish() {
  file_out_->cd();
  config_.Write("Configuration");
  tree_out_->Write();
  file_out_->Close();
  file_in_->Close();
}

void GenericContainerFiller::Run(int nEntries) {
  Init();

  const size_t nTreeEntries = tree_in_->GetEntries();
  const size_t nRunEntries = (nEntries < 0 || nEntries > nTreeEntries) ? nTreeEntries : nEntries;

  int previousTriggerVar{-799};
  for (int iEntry = 0; iEntry < nRunEntries; iEntry++) {
    previousTriggerVar = Exec(iEntry, previousTriggerVar);
  }// iEntry
  tree_out_->Fill();

  Finish();
}

int GenericContainerFiller::DetermineFieldIdByName(const std::vector<IndexMap>& iMap, const std::string& name) {
  auto distance = std::distance(iMap.begin(), std::find_if(iMap.begin(), iMap.end(), [&name](const IndexMap& p) { return p.name_ == name; }));
  if (distance == iMap.size()) throw std::runtime_error("DetermineFieldIdByName(): name " + name + " is missing");
  return distance;
}

void GenericContainerFiller::SetAddressFICS(TBranch* branch, const IndexMap& imap, FICS& ficc) {
  if (imap.field_type_ == "TLeafF") branch->SetAddress(&ficc.float_);
  else if (imap.field_type_ == "TLeafI")
    branch->SetAddress(&ficc.int_);
  else if (imap.field_type_ == "TLeafB")
    branch->SetAddress(&ficc.char_);
  else if (imap.field_type_ == "TLeafS")
    branch->SetAddress(&ficc.short_);
  else
    throw std::runtime_error("GenericContainerFiller::SetAddressFICS(): unsupported filed type " + imap.field_type_);
}

void GenericContainerFiller::SetFieldsFICS(const std::vector<IndexMap>& imap, Container& container, const std::vector<FICS>& ficc) {
  for (int iV = 0; iV < ficc.size(); iV++) {
    if (imap.at(iV).field_type_ == "TLeafF") container.SetField(ficc.at(iV).float_, imap.at(iV).index_);
    else if (imap.at(iV).field_type_ == "TLeafI")
      container.SetField(ficc.at(iV).int_, imap.at(iV).index_);
    else if (imap.at(iV).field_type_ == "TLeafB")
      container.SetField(static_cast<int>(ficc.at(iV).char_), imap.at(iV).index_);
    else if (imap.at(iV).field_type_ == "TLeafS")
      container.SetField(static_cast<int>(ficc.at(iV).short_), imap.at(iV).index_);
    else
      throw std::runtime_error("GenericContainerFiller::SetFieldsFICS(): unsupported filed type " + imap.at(iV).field_type_);
  }
}