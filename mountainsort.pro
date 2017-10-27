TEMPLATE = subdirs

# usage:
# qmake
# qmake "COMPONENTS = mproc"

isEmpty(COMPONENTS) {
    COMPONENTS = mountainsortalg
}

isEmpty(GUI) {
    GUI = on
}

CONFIG += ordered

defineReplace(ifcomponent) {
  contains(COMPONENTS, $$1) {
    message(Enabling $$1)
    return($$2)
  }
  return("")
}

SUBDIRS += $$ifcomponent(mountainsortalg,packages/mountainsortalg/mountainsortalg.pro)

