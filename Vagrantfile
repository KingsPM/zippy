# -*- mode: ruby -*-
# vi: set ft=ruby :
Vagrant.configure(2) do |config|
  #config.vm.synced_folder "..", "/home/vagrant/dev", type: "nfs"  # development folder  
  config.vm.synced_folder "..", "/home/vagrant/dev"
  config.ssh.forward_x11 = true
  config.vm.define "sandbox" do |sandbox|
    sandbox.vm.box = "ubuntu/precise64"
    sandbox.vm.hostname = "sandbox.molpath"
    sandbox.vm.network "private_network", ip: "192.168.33.22"
    sandbox.vm.provider "virtualbox" do |vb|
      vb.gui = false
      vb.cpus = 2
      vb.memory = 2048
    end
    sandbox.vm.provision "shell", path: "minimal.sh"  # base install and dotfile import
    sandbox.vm.provision "docker" do |d|
        d.pull_images "debian:jessie"
    end
  end
end
